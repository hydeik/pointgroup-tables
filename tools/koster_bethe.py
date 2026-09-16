"""Scaffolding for the Koster Bethe-symbol transcription (wayfinder ticket #15).

Three subcommands:

    uv run python tools/koster_bethe.py templates            # write missing template files
    uv run python tools/koster_bethe.py check [GROUP ...]    # validate transcribed files
    uv run python tools/koster_bethe.py generate             # write the shipped Koster table

A template file per crystallographic point group is written to
tests/reference/koster/bethe/ with one entry per irrep (ordinary and
projective) in the canonical irrep order of docs/spec/irrep-labels.md, each
carrying the characters on the identifying elements of the group, so the
transcriber can find the row in Koster's table by its characters and copy
the Bethe symbol without knowing the canonical rules. The checker validates
a filled file: the entries still match the computed irreps, every Bethe
symbol matches the grammar, the parity sign is present exactly in the
centrosymmetric groups and agrees with the g/u suffix, symbols are distinct
within a group, and the numbering is a permutation of 1..n (a warning, since
Koster is the authority).

The prototypes under prototypes/ supply the groups and, for the dihedral
and octahedral families, the irreps; the cyclic and tetrahedral families
are built here from the labelling rules. Nothing here is production code.

``generate`` writes the shipped Koster table
(src/pointgroup_tables/data/conventions/bethe/koster.toml, ticket #30) from
the checked transcriptions: it rebuilds every group with the production
package, verifies that the transcribed symbols and identifying characters
are those the package computes, and records the characters in the value
grammar of the document field so that the package attaches each Bethe
symbol by characters at run time (docs/spec/verification.md, ``REF-8``).
"""

from __future__ import annotations

import argparse
import math
import re
import sys
import tomllib
from fractions import Fraction as Fr
from pathlib import Path

from pointgroup_tables.document import build_model
from pointgroup_tables.field import CyclotomicField, Value
from pointgroup_tables.irreps import build_irreps as production_irreps
from pointgroup_tables.symbols import parse_group
from pointgroup_tables.values import parse as parse_value
from pointgroup_tables.values import pretty as pretty_value

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "prototypes"))
import abstract_group_prototype as ag
import output_schema_prototype as osp

BETHE_DIR = ROOT / "tests" / "reference" / "koster" / "bethe"
OUT_FILE = (
    ROOT
    / "src"
    / "pointgroup_tables"
    / "data"
    / "conventions"
    / "bethe"
    / "koster.toml"
)
SOURCE = (
    "G. F. Koster, J. O. Dimmock, R. G. Wheeler and H. Statz, "
    "Properties of the Thirty-Two Point Groups (Cambridge, MA: MIT Press, 1963)"
)

CRYSTALLOGRAPHIC = [
    "C1", "Ci", "C2", "Cs", "C2h", "D2", "C2v", "D2h",
    "C4", "S4", "C4h", "D4", "C4v", "D2d", "D4h",
    "C3", "S6", "D3", "C3v", "D3d",
    "C6", "C3h", "C6h", "D6", "C6v", "D3h", "D6h",
    "T", "Th", "O", "Td", "Oh",
]  # fmt: skip

BETHE_RE = re.compile(r"^Gamma(\d+)([+-])?$")


# --------------------------------------------------------------------------
# Irreps of the families the output-schema prototype does not build
# --------------------------------------------------------------------------


def _is_centro(G):
    return any(p == -1 and r == G.rotation_group.elements[0] for (r, p) in G.elements)


def _element_by_rotation(G, r, prefer_proper=True):
    cands = [g for g in G.elements if g[0] == r]
    cands.sort(key=lambda g: (g[1] != 1) if prefer_proper else (g[1] == 1))
    return cands[0]


def _zeta_of_order(F, N, p):
    """zeta_N^p as an element of the field (N | M, or N = 2 * odd with odd | M)."""
    M = F.conductor
    if M % N == 0:
        return F.zeta(M // N * p)
    assert N % 2 == 0 and M % (N // 2) == 0, (M, N)
    half = N // 2  # odd: zeta_N = -zeta_half^((half + 1) / 2)
    z = F.zeta(M // half * (p * (half + 1) // 2))
    return z if p % 2 == 0 else F.neg(z)


def _exponent(F, x, N):
    """p in [0, N) with x == zeta_N^p, or None."""
    for p in range(N):
        if x == _zeta_of_order(F, N, p):
            return p
    return None


def cyclic_irreps(ctx):
    """Characters of every irrep of G with R = C_N (Cn, Cnh, S2n, C1, Ci, Cs)."""
    F = ctx.field
    G = ctx.G
    N = int(ctx.R.name[1:])
    signs = (1, -1) if _is_centro(G) else (1,)
    out = []  # (projective, chars)
    # ordinary: chi_j(r) = zeta_N^(j r), times the parity sign in a centro group
    for j in range(N):
        for sign in signs:
            chars = []
            for r, p in G.elements:
                z = _zeta_of_order(F, N, j * r)
                chars.append(F.neg(z) if (p == -1 and sign == -1) else z)
            out.append((False, chars))
    # projective: chi_j(s(r)) = zeta_2N^(j s(r)), j odd, s(r) in Z_2N from the section
    for j in range(1, 2 * N, 2):
        for sign in signs:
            chars = []
            for r, p in G.elements:
                z = _zeta_of_order(F, 2 * N, j * ctx.D.cover.section(r))
                chars.append(F.neg(z) if (p == -1 and sign == -1) else z)
            out.append((True, chars))
    return out


def name_axial(ctx, chars, projective):
    """Canonical Mulliken symbol of an irrep of an axial group from its characters
    (docs/spec/irrep-labels.md, rules 1-11); returns (symbol, partner symbol or None)."""
    F = ctx.field
    G = ctx.G
    idx = {g: k for k, g in enumerate(G.elements)}
    labels = ctx.labels
    dihedral = ctx.R.name.startswith("D")
    N = int(ctx.R.name[1:])
    rot = (lambda k: (k % N, 0)) if dihedral else (lambda k: k % N)
    c = _element_by_rotation(G, rot(1))
    odd_h = ctx.family in ("Ch", "Dh") and ctx.n % 2 == 1
    c0 = _element_by_rotation(G, rot(2)) if odd_h else c
    Nc, Nc0 = ag.element_order(G, c), ag.element_order(G, c0)
    d = int(chars[0][0])
    centro = _is_centro(G)

    def at(label):
        return chars[idx[next(g for g in G.elements if labels[g] == label)]]

    if centro:
        suf = "g" if at("i") == F.rational(d) else "u"
    elif odd_h and not projective:
        suf = "'" if at("sigma_h") == F.rational(d) else "''"
    else:
        suf = ""

    def pair(comp, body):
        other = "2" if comp == "1" else "1"
        return f"{comp}{body}{suf}", f"{other}{body}{suf}"

    if not projective:
        x = chars[idx[c0]]
        if d == 1:
            e = _exponent(F, x, Nc0)
            assert e is not None, (ctx.symbol, x)
            if e != 0 and 2 * e != Nc0:
                k = min(e, Nc0 - e)
                return pair(
                    "1" if e < Nc0 / 2 else "2", "E" + ("" if Nc0 <= 4 else str(k))
                )
            letter = "A" if e == 0 else "B"
            if not dihedral:
                return letter + suf, None
            s = chars[idx[_element_by_rotation(G, (0, 1))]]
            if ctx.n == 2 and ctx.family in ("D", "Dh"):  # rule 6
                xc, xs2 = chars[idx[c]], chars[idx[_element_by_rotation(G, (1, 1))]]
                one = F.rational(1)
                if xc == one and s == one and xs2 == one:
                    return "A" + suf, None
                return ("B1" if xc == one else "B3" if s == one else "B2") + suf, None
            return letter + ("1" if s == F.rational(1) else "2") + suf, None
        for k in range(1, Nc0 // 2 + 1):
            if 2 * k != Nc0 and x == F.add(
                _zeta_of_order(F, Nc0, k), _zeta_of_order(F, Nc0, -k)
            ):
                return "E" + ("" if Nc0 <= 4 else str(k)) + suf, None
        raise AssertionError((ctx.symbol, "ordinary", x))
    x = chars[idx[c]]
    if d == 2:
        for k in range(1, Nc, 2):
            if x == F.add(_zeta_of_order(F, 2 * Nc, k), _zeta_of_order(F, 2 * Nc, -k)):
                return f"E{k}/2{suf}", None
        raise AssertionError((ctx.symbol, "projective", x))
    if dihedral:  # k = N, odd N: the pair with chi(s(s)) = +-i
        s = chars[idx[_element_by_rotation(G, (0, 1))]]
        return pair("1" if s == F.i() else "2", f"E{Nc}/2")
    if x == F.rational(1):
        return "A1/2" + suf, None  # C1, Ci
    if x == F.rational(-1):
        return f"B{Nc}/2{suf}", None  # odd N
    jp = _exponent(F, x, 2 * Nc)
    assert jp is not None and jp % 2 == 1, (ctx.symbol, x, jp)
    k = jp if jp < Nc else 2 * Nc - jp
    return pair("1" if jp < Nc else "2", f"E{k}/2")


def tetrahedral_irreps(ctx):
    """All irreps of T and Th from the quaternion w-coordinate and the C3 classes."""
    F = ctx.field
    G = ctx.G
    centro = _is_centro(G)
    labels = ctx.labels
    zeta3 = _zeta_of_order(F, 3, 1)
    zeta3c = _zeta_of_order(F, 3, 2)
    cls_of = {}
    for cl in ctx.classes:
        for g in cl["members"]:
            cls_of[G.elements[g]] = cl["index"]
    c3 = next(g for g in G.elements if labels[g] == "C3^1(1)")
    c3_class = cls_of[c3]

    def w_of(r):
        a, b = ctx.D.cover.section(r)[0]
        assert b == 0
        return F.rational(a)

    def one_e(r, conjugate):
        proper = (r, 1)
        if ag.element_order(G, proper) != 3:
            return F.rational(1)
        same = cls_of[proper] == c3_class
        return (zeta3c if same else zeta3) if conjugate else (zeta3 if same else zeta3c)

    reps = [
        ("A", False, lambda r: F.rational(1)),
        ("1E", False, lambda r: one_e(r, False)),
        ("2E", False, lambda r: one_e(r, True)),
        (
            "T",
            False,
            lambda r: F.sub(F.scale(F.mul(w_of(r), w_of(r)), Fr(4)), F.rational(1)),
        ),
        ("E1/2", True, lambda r: F.scale(w_of(r), Fr(2))),
        ("1E3/2", True, lambda r: F.mul(F.scale(w_of(r), Fr(2)), one_e(r, False))),
        ("2E3/2", True, lambda r: F.mul(F.scale(w_of(r), Fr(2)), one_e(r, True))),
    ]
    partners = {"1E": "2E", "2E": "1E", "1E3/2": "2E3/2", "2E3/2": "1E3/2"}
    irreps = []
    for symbol, projective, fn in reps:
        for sign in (1, -1) if centro else (1,):
            chars = [
                fn(r) if (p == 1 or sign == 1) else F.neg(fn(r))
                for (r, p) in G.elements
            ]
            suf = ("g" if sign == 1 else "u") if centro else ""
            partner = partners.get(symbol)
            irreps.append(
                osp.make_irrep(
                    ctx,
                    symbol + suf,
                    projective,
                    chars,
                    partner_symbol=(partner + suf) if partner else None,
                )
            )
    return irreps


_original_build = osp.build_irreps


def build_irreps(ctx):
    if ctx.R.name.startswith("C"):
        raw = cyclic_irreps(ctx)
    elif ctx.R.name == "T":
        irreps = tetrahedral_irreps(ctx)
    elif ctx.R.name.startswith("D"):
        raw = [(ir.projective, ir.chars) for ir in _original_build(ctx)]
    else:
        return _original_build(ctx)
    if ctx.R.name != "T":
        irreps = []
        for projective, chars in raw:
            symbol, partner = name_axial(ctx, chars, projective)
            irreps.append(
                osp.make_irrep(ctx, symbol, projective, chars, partner_symbol=partner)
            )
        assert len({ir.symbol for ir in irreps}) == len(irreps), (
            ctx.symbol,
            [ir.symbol for ir in irreps],
        )
    # the tail of the prototype's build_irreps: order, index, partner, Frobenius-Schur
    irreps.sort(key=osp.order_key)
    by_symbol = {}
    for k, ir in enumerate(irreps):
        ir.index = k
        by_symbol[ir.symbol] = k
    F = ctx.field
    for ir in irreps:
        if ir.partner_symbol:
            ir.partner = by_symbol[ir.partner_symbol]
        total = F.zero()
        for a in range(ctx.N):
            sq = ctx.table[a][a]
            term = ir.chars[sq]
            if ir.projective and ctx.omega[a][a] == -1:
                term = F.neg(term)
            total = F.add(total, term)
        total = F.scale(total, Fr(1, ctx.N))
        assert all(x == 0 for x in total[1:]) and total[0].denominator == 1, (
            ir.symbol,
            total,
        )
        ir.frobenius_schur = int(total[0])
    return irreps


# The prototype looks the builder up at call time, so replacing the attribute
# routes every family through the builder above.
osp.build_irreps = build_irreps  # ty: ignore[invalid-assignment]


def context(symbol):
    return osp.Context(symbol)


def self_check(ctx):
    """Orthogonality over G-dagger and the irrep count against the class count."""
    F = ctx.field
    irreps = ctx.irreps
    assert len(irreps) == len(ctx.double_classes), (
        ctx.symbol,
        len(irreps),
        len(ctx.double_classes),
    )
    for a in irreps:
        for b in irreps:
            total = F.zero()
            for k in range(2 * ctx.N):
                x = osp.character_on(ctx, a, k)
                y = F.conj(osp.character_on(ctx, b, k))
                total = F.add(total, F.mul(x, y))
            expect = F.rational(2 * ctx.N if a is b else 0)
            assert total == expect, (ctx.symbol, a.symbol, b.symbol, total)


# --------------------------------------------------------------------------
# Identifying elements and readable values
# --------------------------------------------------------------------------


def identifying_elements(ctx):
    G = ctx.G
    labels = ctx.labels
    out = []

    def add_label(name):
        for g in G.elements:
            if labels[g] == name:
                out.append(g)
                return

    def add_rotation(r):
        out.append(_element_by_rotation(G, r))

    fam, n = ctx.family, ctx.n
    if ctx.polyhedral or ctx.R.name == "T":
        if ctx.R.name == "O":
            add_label("C4^1(1)" if fam == "Oh" or fam == "O" else "S4^-1(1)")
        add_label("C3^1(1)")
        if ctx.R.name == "T":
            add_label("C2(1)")
        add_label("i")
    elif ctx.R.name.startswith("C"):
        N = int(ctx.R.name[1:])
        add_rotation(1 % N)
        if fam == "Ch" and n % 2 == 1:
            add_rotation(2 % N)
        add_label("i")
        add_label("sigma_h")
    else:
        add_rotation((1, 0))
        if fam in ("Dh",) and n % 2 == 1:
            add_rotation((2, 0))
        add_rotation((0, 1))
        if n == 2 and fam in ("D", "Dh"):
            add_rotation((1, 1))
        add_label("i")
        add_label("sigma_h")
    return out


def zeta_name(M, m):
    g = math.gcd(m, M)
    N, k = M // g, m // g
    if N == 1:
        return "1"
    if N == 2:
        return "-1"
    if N == 4:
        return "i" if k == 1 else "-i"
    return f"zeta{N}" if k == 1 else f"zeta{N}^{k}"


def pretty(F, x):
    if all(a == 0 for a in x[1:]):
        return str(x[0])
    M = F.conductor
    found = []
    for m in range(M):
        z = F.zeta(m)
        for r in (1, -1, 2, -2, 3, -3):
            if x == F.scale(z, Fr(r)):
                head = {1: "", -1: "-", 2: "2*", -2: "-2*", 3: "3*", -3: "-3*"}[r]
                found.append(
                    (M // math.gcd(m, M), r < 0, abs(r), head + zeta_name(M, m))
                )
    if found:
        return min(found)[3]
    roots = []
    if M % 8 == 0:
        roots.append(("sqrt2", F.sqrt2()))
    if M % 12 == 0:
        roots.append(("sqrt3", F.add(F.zeta(M // 12), F.zeta(-(M // 12)))))
    if M % 5 == 0:
        roots.append(("phi", F.phi()))
    for name, val in roots:
        for a in range(-3, 4):
            for b in (1, -1, 2, -2):
                if x == F.add(F.rational(a), F.scale(val, Fr(b))):
                    bs = {1: name, -1: f"-{name}", 2: f"2*{name}", -2: f"-2*{name}"}[b]
                    if a == 0:
                        return bs
                    return f"{a} {'+' if b > 0 else '-'} {bs.lstrip('-')}"
    z = F.complex(x)
    return f"{z.real:+.4f}{z.imag:+.4f}i"


def numeric_str(F, x):
    z = F.complex(x)
    re_, im = round(z.real, 4) + 0.0, round(z.imag, 4) + 0.0
    return f"{re_:+.4f}{im:+.4f}i"


def entries(ctx):
    """One (symbol, projective, dimension, [pretty chars], [numeric chars]) per irrep."""
    F = ctx.field
    elems = identifying_elements(ctx)
    idx = {g: k for k, g in enumerate(ctx.G.elements)}
    rows = []
    for ir in ctx.irreps:
        vals = [ir.chars[idx[g]] for g in elems]
        rows.append(
            {
                "symbol": ir.symbol,
                "projective": ir.projective,
                "dimension": ir.dimension,
                "frobenius_schur": ir.frobenius_schur,
                "characters": [pretty(F, v) for v in vals],
                "numeric": [numeric_str(F, v) for v in vals],
                "rational": all(all(a == 0 for a in v[1:]) for v in vals),
            }
        )
    return [ctx.labels[g] for g in elems], rows


# --------------------------------------------------------------------------
# Templates
# --------------------------------------------------------------------------


def toml_str(s):
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'


def toml_list(xs):
    return "[" + ", ".join(toml_str(x) for x in xs) + "]"


def template(symbol):
    ctx = context(symbol)
    self_check(ctx)
    elems, rows = entries(ctx)
    centro = _is_centro(ctx.G)
    n_ord = sum(1 for r in rows if not r["projective"])
    lines = [
        f"# Koster, Dimmock, Wheeler & Statz, Properties of the Thirty-Two Point Groups (1963): Bethe symbols of {symbol}",
        "# Format: tests/reference/koster/bethe/README.md",
        "#",
        f"# {symbol}: {len(rows)} irreps ({n_ord} ordinary, {len(rows) - n_ord} projective), in canonical order.",
        "# Each entry lists the characters on the identifying elements below; for a projective",
        "# irrep they are the characters of the lifted elements in the canonical section.",
        "# Find the row with these characters in Koster's table and copy its Gamma symbol into",
        '# bethe ("Gamma3"'
        + (
            ', with the parity sign in this centrosymmetric group: "Gamma3+" / "Gamma3-").'
            if centro
            else "; no parity sign in this group)."
        ),
        "",
        f"group = {toml_str(symbol)}",
        'table = ""          # Koster\'s table and page for this group, e.g. "Table 33, p. 62"',
        'status = "todo"     # todo | partial | done',
        'notes = ""          # disagreements with Bradley & Cracknell or Bilbao; anything unreadable',
        f"elements = {toml_list(elems)}",
        "",
    ]
    for r in rows:
        lines.append("[[irreps]]")
        lines.append(f"symbol = {toml_str(r['symbol'])}")
        lines.append(f"dimension = {r['dimension']}")
        lines.append(f"projective = {'true' if r['projective'] else 'false'}")
        lines.append(f"characters = {toml_list(r['characters'])}")
        if not r["rational"]:
            lines.append(f"numeric = {toml_list(r['numeric'])}")
        lines.append('bethe = ""')
        lines.append("")
    return "\n".join(lines)


def write_templates(force=False):
    BETHE_DIR.mkdir(parents=True, exist_ok=True)
    written = []
    for symbol in CRYSTALLOGRAPHIC:
        path = BETHE_DIR / f"{ag.ALIASES.get(symbol, symbol)}.toml"
        if path.exists() and not force:
            continue
        path.write_text(template(symbol) + "\n")
        written.append(path.name)
    print(f"wrote {len(written)} template(s): {' '.join(written) or '(none missing)'}")


# --------------------------------------------------------------------------
# Checker
# --------------------------------------------------------------------------


class Report:
    def __init__(self, name):
        self.name = name
        self.errors = []
        self.warnings = []

    def error(self, msg):
        self.errors.append(msg)

    def warn(self, msg):
        self.warnings.append(msg)


def check_file(path):
    rep = Report(path.name)
    data = tomllib.loads(path.read_text())
    symbol = data.get("group")
    if symbol not in CRYSTALLOGRAPHIC and symbol not in ag.ALIASES:
        rep.error(f"group {symbol!r} is not one of the 32 crystallographic groups")
        return rep
    ctx = context(symbol)
    elems, rows = entries(ctx)
    centro = _is_centro(ctx.G)
    if data.get("elements") != elems:
        rep.error(
            f"elements {data.get('elements')} differ from the computed {elems}; regenerate"
        )
    file_rows = data.get("irreps", [])
    if [r.get("symbol") for r in file_rows] != [r["symbol"] for r in rows]:
        rep.error(
            "irrep symbols or order differ from the computed ones; regenerate the template"
        )
        return rep
    for fr, cr in zip(file_rows, rows):
        if fr.get("characters") != cr["characters"]:
            rep.error(
                f"{cr['symbol']}: characters {fr.get('characters')} differ from computed {cr['characters']}"
            )
    filled = [(cr, fr.get("bethe", "")) for fr, cr in zip(file_rows, rows)]
    n_filled = sum(1 for _, b in filled if b)
    status = data.get("status")
    expected = (
        "todo" if n_filled == 0 else ("done" if n_filled == len(rows) else "partial")
    )
    if status != expected:
        rep.warn(
            f"status is {status!r} but {n_filled}/{len(rows)} entries are filled ({expected!r})"
        )
    seen = {}
    numbers = {"+": [], "-": [], "": []}
    by_symbol = {cr["symbol"]: b for cr, b in filled}
    for cr, b in filled:
        if not b:
            continue
        m = BETHE_RE.match(b)
        if not m:
            rep.error(f"{cr['symbol']}: bethe {b!r} does not match Gamma<n>[+-]")
            continue
        num, sign = int(m.group(1)), m.group(2) or ""
        if centro and not sign:
            rep.error(
                f"{cr['symbol']}: {b} needs a parity sign in a centrosymmetric group"
            )
        if not centro and sign:
            rep.error(
                f"{cr['symbol']}: {b} carries a parity sign but {symbol} is not centrosymmetric"
            )
        if centro and sign:
            suf = "g" if cr["symbol"].endswith("g") else "u"
            if (sign == "+") != (suf == "g"):
                rep.error(
                    f"{cr['symbol']}: sign of {b} disagrees with the {suf} suffix"
                )
        if b in seen:
            rep.error(f"{cr['symbol']}: {b} already used by {seen[b]}")
        seen[b] = cr["symbol"]
        numbers[sign].append((num, cr))
        partner = None
        if cr["symbol"][0] in "12":
            partner = ("2" if cr["symbol"][0] == "1" else "1") + cr["symbol"][1:]
        if partner and by_symbol.get(partner) == b:
            rep.error(f"{cr['symbol']} and its partner {partner} share {b}")
    if n_filled == len(rows):
        for sign, lst in numbers.items():
            if not lst:
                continue
            nums = sorted(n for n, _ in lst)
            if nums != list(range(1, len(nums) + 1)):
                rep.warn(
                    f"numbers{sign or ''} {nums} are not 1..{len(nums)}; fine if Koster does that, note it"
                )
            ord_max = max((n for n, cr in lst if not cr["projective"]), default=0)
            proj_min = min((n for n, cr in lst if cr["projective"]), default=10**6)
            if proj_min < ord_max:
                rep.warn(
                    f"a projective irrep has a smaller number{sign} than an ordinary one"
                )
        if centro:
            plus = {n: cr["symbol"] for n, cr in numbers["+"]}
            for n, cr in numbers["-"]:
                s = plus.get(n)
                if s is not None and s[:-1] != cr["symbol"][:-1]:
                    rep.warn(
                        f"Gamma{n}+ is {s} but Gamma{n}- is {cr['symbol']}; expected the g/u pair"
                    )
    return rep


def run_check(groups):
    paths = sorted(BETHE_DIR.glob("*.toml"))
    if groups:
        wanted = {ag.ALIASES.get(g, g) for g in groups}
        paths = [p for p in paths if p.stem in wanted]
    failed = 0
    for path in paths:
        rep = check_file(path)
        data = tomllib.loads(path.read_text())
        print(
            f"{path.name}: status={data.get('status')} errors={len(rep.errors)} warnings={len(rep.warnings)}"
        )
        for e in rep.errors:
            print(f"  error: {e}")
        for w in rep.warnings:
            print(f"  warning: {w}")
        failed += bool(rep.errors)
    print(f"{len(paths)} file(s), {failed} with errors")
    return 1 if failed else 0


# --------------------------------------------------------------------------
# The shipped Koster table
# --------------------------------------------------------------------------

_TRANSCRIBED_ROOT = re.compile(r"^(-)?(?:(\d+)\*)?zeta(\d+)(?:\^(\d+))?$")


def transcribed_value(F: CyclotomicField, text: str) -> Value:
    """A character as the transcription spells it: the template's ``zetaN^k``
    (a root of unity named by its own order) or a string of the value grammar."""
    match = _TRANSCRIBED_ROOT.match(text)
    if match is None:
        return parse_value(F, text)
    sign, factor, order, power = match.groups()
    value = F.root_of_unity(int(order), int(power or 1))
    if factor:
        value = F.scale(value, Fr(int(factor)))
    return F.neg(value) if sign else value


class GenerateError(Exception):
    """A transcription the production package does not reproduce."""


def shipped_group(symbol: str) -> tuple[str, dict]:
    """The shipped entry of one group, verified against the production irreps."""
    spec = parse_group(symbol)
    path = BETHE_DIR / f"{spec.symbol}.toml"
    data = tomllib.loads(path.read_text())
    if data.get("status") != "done":
        raise GenerateError(f"{path.name}: status is {data.get('status')!r}, not done")
    model = build_model(spec)
    F = model.field
    by_label = {label.text: k for k, label in enumerate(model.labels)}
    missing = [name for name in data["elements"] if name not in by_label]
    if missing:
        raise GenerateError(f"{path.name}: elements {missing} are not canonical labels")
    positions = [by_label[name] for name in data["elements"]]
    # the raw irreps, not ``model.irreps``: that property attaches the very
    # table this command writes
    computed = production_irreps(
        model.group, model.double, F, model.table, model.factor
    )
    symbols = [ir.symbol.text for ir in computed]
    if [row["symbol"] for row in data["irreps"]] != symbols:
        raise GenerateError(
            f"{path.name}: irrep symbols or order differ from the package's {symbols}"
        )
    rows = []
    keys = set()
    for ir, row in zip(computed, data["irreps"], strict=True):
        values = [ir.characters[k] for k in positions]
        expected = [transcribed_value(F, text) for text in row["characters"]]
        if values != expected:
            raise GenerateError(
                f"{path.name}: {ir.symbol.text} has characters "
                f"{[pretty_value(F, v) for v in values]}, the transcription says {row['characters']}"
            )
        if ir.projective != row["projective"]:
            raise GenerateError(
                f"{path.name}: {ir.symbol.text} projective flag differs"
            )
        characters = [pretty_value(F, v) for v in values]
        key = (ir.projective, tuple(characters))
        if key in keys:
            raise GenerateError(
                f"{path.name}: the identifying elements do not separate {ir.symbol.text}"
            )
        keys.add(key)
        rows.append(
            {
                "symbol": ir.symbol.text,
                "projective": ir.projective,
                "characters": characters,
                "bethe": row["bethe"],
            }
        )
    entry = {
        "table": data["table"],
        "conductor": F.conductor,
        "elements": list(data["elements"]),
        "irreps": rows,
    }
    return spec.symbol, entry


def emit_shipped(groups: list[tuple[str, dict]]) -> str:
    lines = [
        "# The Koster table: Bethe symbols of the irreps of the 32 crystallographic point",
        '# groups (docs/spec/irrep-labels.md, "Bethe aliases"). Generated by',
        "# `uv run python tools/koster_bethe.py generate` from the transcriptions under",
        "# tests/reference/koster/bethe/; do not edit by hand.",
        "#",
        "# Each group lists its identifying elements (canonical labels) and, per irrep, the",
        "# characters on them in the value grammar of the group's canonical field (for a",
        "# projective irrep, of the lifted elements in the canonical section). The package",
        "# attaches a Bethe symbol to the computed irrep with those characters (REF-8).",
        "",
        f"name = {toml_str(KOSTER_NAME)}",
        f"source = {toml_str(SOURCE)}",
        "",
    ]
    for symbol, entry in groups:
        lines += [
            f"[groups.{symbol}]",
            f"table = {toml_str(entry['table'])}",
            f"conductor = {entry['conductor']}",
            f"elements = {toml_list(entry['elements'])}",
            "",
        ]
        for row in entry["irreps"]:
            lines += [
                f"[[groups.{symbol}.irreps]]",
                f"symbol = {toml_str(row['symbol'])}",
                f"projective = {'true' if row['projective'] else 'false'}",
                f"characters = {toml_list(row['characters'])}",
                f"bethe = {toml_str(row['bethe'])}",
                "",
            ]
    return "\n".join(lines)


KOSTER_NAME = "koster"


def generate() -> int:
    """Check every transcription, then write the shipped Koster table."""
    if run_check([]) != 0:
        print("not generated: the transcriptions do not pass the checker")
        return 1
    groups = []
    try:
        for symbol in CRYSTALLOGRAPHIC:
            groups.append(shipped_group(symbol))
    except GenerateError as exc:
        print(f"not generated: {exc}")
        return 1
    OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    OUT_FILE.write_text(emit_shipped(groups) + "\n")
    print("wrote", OUT_FILE.relative_to(ROOT), f"({len(groups)} groups)")
    return 0


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    t = sub.add_parser("templates", help="write missing template files")
    t.add_argument("--force", action="store_true", help="overwrite existing files")
    c = sub.add_parser("check", help="validate transcribed files")
    c.add_argument("groups", nargs="*")
    s = sub.add_parser("show", help="print the computed entries of a group")
    s.add_argument("group")
    sub.add_parser(
        "generate", help="write the shipped Koster table from the checked files"
    )
    args = ap.parse_args(argv)
    if args.cmd == "templates":
        write_templates(force=args.force)
        return 0
    if args.cmd == "generate":
        return generate()
    if args.cmd == "show":
        print(template(args.group))
        return 0
    return run_check(args.groups)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
