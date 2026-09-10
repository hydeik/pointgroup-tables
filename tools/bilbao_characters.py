"""Structural comparison of the Bilbao character tables with the computed
tables (the prototype of ``REF-2``, wayfinder ticket #20).

    uv run python tools/bilbao_characters.py check [GROUP ...]

For every group the checker reads ``tests/reference/bilbao/characters/
<group>.toml`` (the POINT page: classes of G, ordinary irreps) and
``<group>-double.toml`` (the REPRESENTATIONS DPG page: classes of G†, every
irrep), rebuilds the group with the prototypes, and searches for a bijection
of classes (candidates by size, element order, parity and, for G†, whether
the class holds lifted, barred or both kinds of elements) under which every
row of the server's table equals one computed irrep exactly. The double-group
file needs no sign freedom: the server lists the classes of G† itself. The
server's reality (1 real, 0 complex, -1 pseudoreal) is compared with the
Frobenius-Schur indicator. After the match the checker reports how the
server's labels sit on the canonical ones: which member of a complex pair
the server calls ``1E``, and where the ``Gamma`` numbers of the two pages
disagree with each other (the ``C6`` case). Nothing here is production code.
"""

from __future__ import annotations

import itertools
import math
import re
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "prototypes"))
sys.path.insert(0, str(ROOT / "tools"))
import abstract_group_prototype as ag
import canonical_basis_prototype as cb
import koster_bethe  # noqa: F401  (extends the schema prototype to the cyclic and T families)
import output_schema_prototype as osp

CHAR_DIR = ROOT / "tests" / "reference" / "bilbao" / "characters"
FS_OF_REALITY = {1: 1, 0: 0, -1: -1}


def table_order(table, identity, g):
    n, h = 1, g
    while h != identity:
        h = table[h][g]
        n += 1
    return n


def identity_of(table):
    return next(
        e for e in range(len(table)) if all(table[e][x] == x for x in range(len(table)))
    )


class Computed:
    """Classes of G and G† with their invariants, and the irreps."""

    def __init__(self, symbol):
        self.ctx = ctx = osp.Context(symbol)
        self.F = ctx.field
        N = ctx.N
        e = identity_of(ctx.table)
        self.single = []  # (invariants, representative, label)
        for cl in ctx.classes:
            rep = cl["representative"]
            inv = (cl["size"], table_order(ctx.table, e, rep), ctx.G.elements[rep][1])
            self.single.append((inv, rep, cl["label"]))
        e2 = identity_of(ctx.double_table)
        self.double = []
        for cl in ctx.double_classes:
            rep = cl["representative"]
            kind = "mixed" if not cl["split"] else ("all" if cl["barred"] else "none")
            inv = (
                cl["size"],
                table_order(ctx.double_table, e2, rep),
                ctx.G.elements[rep % N][1],
                kind,
            )
            self.double.append((inv, rep, cl["label"]))
        self.irreps = ctx.irreps

    def char(self, ir, q):
        """Character of irrep ``ir`` on the double-group element ``q``."""
        N = self.ctx.N
        x = ir.chars[q % N]
        return self.F.neg(x) if (q >= N and ir.projective) else x


def lift_to(F, G, x):
    return x if F.conductor == G.conductor else cb.lift(F, G, x)


SEITZ_RE = re.compile(
    r"^(?P<d>d)?(?P<imp>-)?(?P<body>[12346]|m)(?P<sense>[+-])?(?:_(?P<dir>.*))?$"
)


def preferred(book_label, comp_label):
    """Geometric tie-break between classes of equal invariants: the server's
    ``n+`` (about 001, or ``C<k>`` holding it) is the rotation by +2 pi / n,
    the canonical ``Cn``; its rotoinversion ``-n+`` is ``i Cn``, the canonical
    ``S^-1`` (rule 3 of the label conventions)."""
    m = SEITZ_RE.match(book_label.split()[-1])  # DPG labels carry their first element
    if not m or not m.group("sense"):
        return False
    positive = (m.group("sense") == "+") != bool(m.group("imp"))
    return positive == ("^-" not in comp_label)


def structural_match(book_classes, book_rows, comp_classes, comp_rows, kinds=None):
    """book_classes: [(invariants, label)]; book_rows: [(label, values, projective)];
    comp_classes: [(invariants, label)]; comp_rows: [(symbol, values, projective)],
    values in one field. ``kinds`` (double group only) gives the kind of
    operations of each book and computed class ("none", "all", "mixed"
    barred), a tie-break only: an isomorphism of G† may exchange a lifted
    class with its barred partner, and the server's section does so on some
    binary rotations and planes. Returns (perm, pairing) or None."""
    n = len(book_classes)
    if sorted(i for i, _ in book_classes) != sorted(i for i, _ in comp_classes):
        return None
    book_kind, comp_kind = kinds or ([None] * n, [None] * n)
    cands = [
        sorted(
            (j for j in range(n) if comp_classes[j][0] == book_classes[i][0]),
            key=lambda j: (
                not preferred(book_classes[i][1], comp_classes[j][1]),
                book_kind[i] != comp_kind[j],
            ),
        )
        for i in range(n)
    ]
    order_ = sorted(range(n), key=lambda i: len(cands[i]))
    proj_comp = [proj for _, _, proj in comp_rows]

    def feasible(perm, assigned):
        """Every book row still has a computed row agreeing on the assigned classes."""
        for _, vals, proj in book_rows:
            if not any(
                proj_comp[s] == proj and all(row[perm[i]] == vals[i] for i in assigned)
                for s, (_, row, _) in enumerate(comp_rows)
            ):
                return False
        return True

    def try_rows(perm):
        used, pairing = set(), []
        for _, vals, proj in book_rows:
            target = None
            for s, (_, row, cproj) in enumerate(comp_rows):
                if s in used or cproj != proj:
                    continue
                if all(row[perm[i]] == vals[i] for i in range(n)):
                    target = s
                    break
            if target is None:
                return None
            used.add(target)
            pairing.append(target)
        return pairing

    def search(pos, perm, used):
        if pos == n:
            pairing = try_rows(perm)
            return (perm[:], pairing) if pairing is not None else None
        i = order_[pos]
        for j in cands[i]:
            if j in used:
                continue
            perm[i] = j
            used.add(j)
            if feasible(perm, order_[: pos + 1]):
                found = search(pos + 1, perm, used)
                if found:
                    return found
            used.discard(j)
        return None

    return search(0, [None] * n, set())


def read_file(path):
    data = tomllib.loads(path.read_text())
    F = cb.field(data["conductor"])
    rows = []
    for r in data["irreps"]:
        rows.append(
            (
                r["label"],
                [cb.parse_value(F, s) for s in r["characters"]],
                r["projective"],
                r,
            )
        )
    return data, F, rows


def check_file(path, comp, double):
    data, Fb, rows = read_file(path)
    F = cb.field(math.lcm(Fb.conductor, comp.F.conductor))
    kinds = None
    if double:
        book_classes = [
            ((c["size"], c["order"], c["parity"]), c["label"] + " " + c["elements"][0])
            for c in data["classes"]
        ]
        comp_classes = [(inv[:3], label) for inv, _, label in comp.double]
        kinds = (
            [c["barred"] for c in data["classes"]],
            [inv[3] for inv, _, _ in comp.double],
        )
        reps = [rep for _, rep, _ in comp.double]
        comp_rows = [
            (
                ir.symbol,
                [lift_to(comp.F, F, comp.char(ir, q)) for q in reps],
                ir.projective,
            )
            for ir in comp.irreps
        ]
    else:
        book_classes = [
            ((c["size"], c["order"], c["parity"]), c["label"]) for c in data["classes"]
        ]
        comp_classes = [(inv, label) for inv, _, label in comp.single]
        reps = [rep for _, rep, _ in comp.single]
        comp_rows = [
            (ir.symbol, [lift_to(comp.F, F, ir.chars[g]) for g in reps], ir.projective)
            for ir in comp.irreps
            if not ir.projective
        ]
    book_rows = [
        (lab, [lift_to(Fb, F, v) for v in vals], proj) for lab, vals, proj, _ in rows
    ]
    lines = []
    ok = True
    if sorted(i for i, _ in book_classes) != sorted(i for i, _ in comp_classes):
        lines.append(
            f"  class invariants differ: server {sorted(i for i, _ in book_classes)} vs computed {sorted(i for i, _ in comp_classes)}"
        )
        return lines, False, None
    found = structural_match(book_classes, book_rows, comp_classes, comp_rows, kinds)
    if not found:
        lines.append(
            "  NO structural match: no class/irrep bijection reproduces the server's table"
        )
        return lines, False, None
    perm, pairing = found
    lines.append(
        "  classes: "
        + ", ".join(
            f"{book_classes[i][1]} = {comp_classes[perm[i]][1]}"
            for i in range(len(perm))
        )
    )
    if kinds:
        swapped = [
            f"{book_classes[i][1]} = {comp_classes[perm[i]][1]}"
            for i in range(len(perm))
            if kinds[0][i] != kinds[1][perm[i]]
        ]
        if swapped:
            lines.append(
                "  section: the server's lifted operations are canonical barred ones (and vice versa) in "
                + ", ".join(swapped)
            )
    lines.append(
        "  irreps:  "
        + ", ".join(
            f"{book_rows[r][0]} = {comp_rows[pairing[r]][0]}"
            for r in range(len(pairing))
        )
    )
    if len(pairing) != len(comp_rows):
        lines.append(
            f"  the server lists {len(pairing)} irreps, the computed table has {len(comp_rows)}"
        )
        ok = False
    symbol_of = {}
    for r, (lab, _, _, rec) in enumerate(rows):
        sym = comp_rows[pairing[r]][0]
        symbol_of[sym] = rec
        if double:
            ir = next(ir for ir in comp.irreps if ir.symbol == sym)
            if FS_OF_REALITY[rec["reality"]] != ir.frobenius_schur:
                lines.append(
                    f"  reality of {lab} is {rec['reality']} but the Frobenius-Schur indicator of {sym} is {ir.frobenius_schur}"
                )
                ok = False
    return lines, ok, symbol_of


def label_report(group, single, double):
    """How the server's labels sit on the canonical symbols."""
    lines = []
    for name, page in (("POINT", single), ("DPG", double)):
        swapped = sorted(
            f"{rec['label']} = {sym}"
            for sym, rec in page.items()
            if sym[0] in "12" and rec["label"][0] in "12" and rec["label"][0] != sym[0]
        )
        same = sorted(
            rec["label"]
            for sym, rec in page.items()
            if sym[0] in "12" and rec["label"][0] == sym[0]
        )
        if swapped:
            lines.append(
                f"  1E/2E on the {name} page: opposite exponent sign to the canonical symbol for {', '.join(swapped)}"
            )
        if same:
            lines.append(
                f"  1E/2E on the {name} page: same exponent sign as the canonical symbol for {', '.join(same)}"
            )
    return lines


AXES = {"x": "100", "y": "010", "z": "001"}


def class_candidates(point_label, dpg_classes):
    """Double-group classes that hold the operation of a POINT class header
    (same type, sense and direction; a header without direction prefers the
    principal axis 001; a direction absent from the page, as for the two
    settings of C3v, is ignored)."""
    m = SEITZ_RE.match(point_label)
    assert m, point_label
    want_dir = AXES.get(m.group("dir") or "", m.group("dir"))

    def hits(strict):
        out = []
        for k, c in enumerate(dpg_classes):
            for e in c["elements"]:
                em = SEITZ_RE.match(e)
                assert em, e
                if (
                    em.group("d")
                    or em.group("body") != m.group("body")
                    or bool(em.group("imp")) != bool(m.group("imp"))
                ):
                    continue
                if m.group("sense") and em.group("sense") != m.group("sense"):
                    continue
                if strict and want_dir and em.group("dir") != want_dir:
                    continue
                out.append((k, em.group("dir")))
                break
        return out

    found = hits(True) or hits(False)
    found.sort(key=lambda kd: kd[1] != "001")
    return [k for k, _ in found]


def page_consistency(single_data, double_data):
    """Match the POINT classes to the DPG classes by their operations and the
    ordinary irreps by their characters; report label and Gamma disagreements."""
    Fs, Fd = cb.field(single_data["conductor"]), cb.field(double_data["conductor"])
    F = cb.field(math.lcm(Fs.conductor, Fd.conductor))
    cands = [
        class_candidates(c["label"], double_data["classes"])
        for c in single_data["classes"]
    ]
    if any(not c for c in cands):
        return ["  page comparison impossible: a POINT class has no DPG class"]
    srows = [
        (r, [lift_to(Fs, F, cb.parse_value(Fs, v)) for v in r["characters"]])
        for r in single_data["irreps"]
    ]
    drows = [
        (r, [lift_to(Fd, F, cb.parse_value(Fd, v)) for v in r["characters"]])
        for r in double_data["irreps"]
        if not r["projective"]
    ]
    for choice in itertools.product(*cands):
        if len(set(choice)) != len(choice):
            continue
        pairs = []
        for srec, svals in srows:
            hit = [
                drec
                for drec, dvals in drows
                if all(dvals[k] == v for k, v in zip(choice, svals))
            ]
            if len(hit) != 1:
                break
            pairs.append((srec, hit[0]))
        else:
            lines = [
                "  POINT classes on the DPG page: "
                + ", ".join(
                    f"{c['label']} = {double_data['classes'][k]['label']}"
                    for c, k in zip(single_data["classes"], choice)
                )
            ]
            diff = [
                f"POINT {s['label']} {s['bethe']} is DPG {d['label']} {d['bethe']}"
                for s, d in pairs
                if s["label"] != d["label"] or s["bethe"] != d["bethe"]
            ]
            if diff:
                lines.append(
                    "  the two pages disagree on ordinary irreps: " + "; ".join(diff)
                )
            else:
                lines.append(
                    "  both pages give the same Mulliken and Gamma labels to the ordinary irreps"
                )
            return lines
    return [
        "  page comparison impossible: no class matching makes every POINT row a DPG row"
    ]


def check(groups):
    wanted = [ag.ALIASES.get(g, g) for g in groups]
    paths = sorted(p for p in CHAR_DIR.glob("*.toml") if not p.stem.endswith("-double"))
    if wanted:
        paths = [p for p in paths if p.stem in wanted]
    bad = 0
    for path in paths:
        symbol = path.stem
        comp = Computed(symbol)
        lines1, ok1, single = check_file(path, comp, double=False)
        lines2, ok2, double = check_file(
            path.with_name(f"{symbol}-double.toml"), comp, double=True
        )
        ok = ok1 and ok2
        print(f"{symbol}: {'ok' if ok else 'MISMATCH'}")
        print(" POINT page")
        for ln in lines1:
            print(ln)
        print(" REPRESENTATIONS DPG page")
        for ln in lines2:
            print(ln)
        if single is not None and double is not None:
            for ln in label_report(symbol, single, double):
                print(ln)
            single_data = tomllib.loads(path.read_text())
            double_data = tomllib.loads(
                path.with_name(f"{symbol}-double.toml").read_text()
            )
            for ln in page_consistency(single_data, double_data):
                print(ln)
        bad += not ok
    print(f"{len(paths)} group(s), {bad} with mismatches")
    return 1 if bad else 0


def main(argv):
    if argv and argv[0] == "check":
        return check(argv[1:])
    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
