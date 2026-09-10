"""Prototype: output document layout of every table (wayfinder ticket #11).

Throwaway code that produces *sample documents* in the layout proposed for
docs/spec/output-schema.md, so the layout can be judged on real data:

* ``elements``, ``cayley`` and ``factor`` for any group the two earlier
  prototypes build (canonical order, canonical labels, conjugacy classes of
  G and G-dagger, the spin factor system);
* ``characters`` for C3v, D6h (axial: all irreps from the canonical-basis
  formulas) and Oh, Ih (polyhedral: characters from the quaternion ``w``
  coordinate of the canonical section, no matrices yet -- those are ticket
  #16);
* ``matrices`` for the axial groups only.

Values are exact in one cyclotomic field per document (docs/spec/values.md);
numeric approximations are plain ``complex`` arithmetic here, not the
mpmath-rounded values of the spec. Bethe symbols are omitted (ticket #15).

    uv run python prototypes/output_schema_prototype.py selfcheck
    uv run python prototypes/output_schema_prototype.py characters C3v
    uv run python prototypes/output_schema_prototype.py all D6h --irreps projective --values exact
    uv run python prototypes/output_schema_prototype.py elements Oh --double-group -o oh.json
    uv run python prototypes/output_schema_prototype.py samples   # writes prototypes/samples/

Nothing here is production code.
"""

from __future__ import annotations

import argparse
import cmath
import json
import math
import re
import sys
from fractions import Fraction as Fr
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import abstract_group_prototype as ag
import double_group_prototype as dg

PROTOTYPE_VERSION = "0.0.0-prototype"

# --------------------------------------------------------------------------
# Exact cyclotomic field Q(zeta_M), elements as coefficient tuples over the
# power basis 1, zeta, ..., zeta^(phi(M)-1)   (docs/spec/values.md)
# --------------------------------------------------------------------------


def _divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]


def _poly_div_exact(num, den):
    """Exact division of integer polynomials (coefficient lists, low to high)."""
    num = list(num)
    q = [0] * (len(num) - len(den) + 1)
    for k in range(len(num) - len(den), -1, -1):
        c = num[k + len(den) - 1] // den[-1]
        q[k] = c
        for j, d in enumerate(den):
            num[k + j] -= c * d
    assert all(c == 0 for c in num), "not an exact division"
    return q


_CYCLO = {}


def cyclotomic_polynomial(M):
    if M not in _CYCLO:
        poly = [-1] + [0] * (M - 1) + [1]  # x^M - 1
        for d in _divisors(M):
            if d < M:
                poly = _poly_div_exact(poly, cyclotomic_polynomial(d))
        _CYCLO[M] = poly
    return _CYCLO[M]


class Field:
    def __init__(self, M):
        self.conductor = M
        self.poly = cyclotomic_polynomial(M)
        self.degree = len(self.poly) - 1

    # construction
    def reduce(self, coeffs):
        c = [Fr(x) for x in coeffs]
        n = self.degree
        for k in range(len(c) - 1, n - 1, -1):
            if c[k]:
                q = c[k]
                for j, p in enumerate(self.poly):
                    c[k - n + j] -= q * p
        c = c[:n] + [Fr(0)] * (n - len(c))
        return tuple(c[:n])

    def rational(self, r):
        return self.reduce([Fr(r)])

    def zeta(self, k=1):
        k %= self.conductor
        return self.reduce([Fr(0)] * k + [Fr(1)])

    def i(self):
        assert self.conductor % 4 == 0
        return self.zeta(self.conductor // 4)

    def sqrt2(self):
        assert self.conductor % 8 == 0
        z = self.conductor // 8
        return self.add(self.zeta(z), self.zeta(-z))

    def sqrt5(self):
        assert self.conductor % 5 == 0
        z = self.conductor // 5
        return self.add(
            self.sub(self.zeta(z), self.zeta(2 * z)),
            self.sub(self.zeta(4 * z), self.zeta(3 * z)),
        )

    def phi(self):
        return self.scale(self.add(self.rational(1), self.sqrt5()), Fr(1, 2))

    # arithmetic
    def add(self, x, y):
        return tuple(a + b for a, b in zip(x, y))

    def sub(self, x, y):
        return tuple(a - b for a, b in zip(x, y))

    def neg(self, x):
        return tuple(-a for a in x)

    def scale(self, x, r):
        return tuple(a * r for a in x)

    def mul(self, x, y):
        n = self.degree
        prod = [Fr(0)] * (2 * n - 1)
        for i, a in enumerate(x):
            if a:
                for j, b in enumerate(y):
                    prod[i + j] += a * b
        return self.reduce(prod)

    def power(self, x, k):
        out = self.rational(1)
        for _ in range(k):
            out = self.mul(out, x)
        return out

    def conj(self, x):
        out = self.rational(0)
        for j, a in enumerate(x):
            if a:
                out = self.add(out, self.scale(self.zeta(-j), a))
        return out

    def zero(self):
        return self.rational(0)

    def is_zero(self, x):
        return all(a == 0 for a in x)

    # output forms
    def exact(self, x):
        return [[a.numerator, a.denominator] for a in x]

    def complex(self, x):
        z = 0j
        for j, a in enumerate(x):
            if a:
                z += float(a) * cmath.exp(2j * math.pi * j / self.conductor)
        return z

    def numeric(self, x):
        z = self.complex(x)
        re_, im = float(z.real), float(z.imag)
        return [0.0 if re_ == 0 else re_, 0.0 if im == 0 else im]

    def metadata(self):
        return {
            "kind": "cyclotomic",
            "conductor": self.conductor,
            "degree": self.degree,
            "generator": f"exp(2*pi*i/{self.conductor})",
            "basis": "power",
        }


# matrices over a field: lists of rows
def mat_mul(F, A, B):
    n, m, p = len(A), len(B), len(B[0])
    return [
        [_sum(F, (F.mul(A[i][k], B[k][j]) for k in range(m))) for j in range(p)]
        for i in range(n)
    ]


def _sum(F, terms):
    out = F.zero()
    for t in terms:
        out = F.add(out, t)
    return out


def mat_scale(F, A, s):
    return [[F.mul(a, s) for a in row] for row in A]


def mat_identity(F, d):
    return [[F.rational(1 if i == j else 0) for j in range(d)] for i in range(d)]


def mat_power(F, A, k):
    out = mat_identity(F, len(A))
    for _ in range(k):
        out = mat_mul(F, out, A)
    return out


def mat_trace(F, A):
    return _sum(F, (A[i][i] for i in range(len(A))))


def mat_exact(F, A):
    return [[F.exact(a) for a in row] for row in A]


def mat_numeric(F, A):
    return [[F.numeric(a) for a in row] for row in A]


# --------------------------------------------------------------------------
# Structured element records from canonical labels
# (docs/spec/label-conventions.md, "Structured element record")
# --------------------------------------------------------------------------

_LABEL = re.compile(
    r"^(?P<bar>~)?(?:(?P<E>E)|(?P<i>i)"
    r"|(?P<C2>C2)(?:\((?P<C2k>\d+)\))?"
    r"|C2(?P<primes>''?)\((?P<seck>\d+)\)"
    r"|(?P<CS>[CS])(?P<m>\d+)\^(?P<j>-?\d+)(?:\((?P<k>\d+)\))?"
    r"|sigma_h(?:\((?P<hk>\d+)\))?"
    r"|sigma_(?P<vd>[vd])\((?P<vdk>\d+)\)"
    r"|sigma\((?P<pk>\d+)\))$"
)


def element_record(index, label, parity, polyhedral):
    m = _LABEL.match(label)
    assert m, label
    g = m.groupdict()
    rec = {"index": index, "label": label, "parity": parity}
    axis_class = axis_index = None
    if g["E"]:
        rec.update(kind="identity", order=1)
    elif g["i"]:
        rec.update(kind="inversion", order=2)
    elif g["C2"]:
        rec.update(kind="rotation", order=2)
        axis_class, axis_index = (
            ("polyhedral", g["C2k"]) if g["C2k"] else ("principal", None)
        )
    elif g["primes"] is not None:
        rec.update(kind="rotation", order=2)
        axis_class = "secondary" if g["primes"] == "'" else "secondary_double_prime"
        axis_index = g["seck"]
    elif g["CS"]:
        rec.update(
            kind="rotation" if g["CS"] == "C" else "rotoreflection", order=int(g["m"])
        )
        rec["power"] = int(g["j"])
        axis_class, axis_index = (
            ("polyhedral", g["k"]) if g["k"] else ("principal", None)
        )
    elif label.startswith("sigma_h"):
        rec.update(kind="reflection", order=2)
        axis_class, axis_index = (
            ("polyhedral", g["hk"]) if g["hk"] else ("principal", None)
        )
    elif g["vd"]:
        rec.update(kind="reflection", order=2)
        # illustrative: the plane is indexed by the secondary axis it contains
        axis_class = "secondary" if g["vd"] == "v" else "secondary_double_prime"
        axis_index = g["vdk"]
    else:
        rec.update(kind="reflection", order=2)
        axis_class, axis_index = "polyhedral", g["pk"]
    if axis_class:
        rec["axis_class"] = axis_class
    if axis_index is not None:
        rec["axis_index"] = int(axis_index)
    return rec


# --------------------------------------------------------------------------
# Conjugacy classes and class labels (docs/spec/conjugacy-classes.md)
# --------------------------------------------------------------------------


def inverses(table):
    n = len(table)
    return [next(j for j in range(n) if table[i][j] == 0) for i in range(n)]


def classes_from_table(table):
    """Classes of a group given by its Cayley table, ordered by first member."""
    n = len(table)
    inv = inverses(table)
    seen = [False] * n
    out = []
    for g in range(n):
        if seen[g]:
            continue
        cls = sorted({table[table[h][g]][inv[h]] for h in range(n)})
        for x in cls:
            seen[x] = True
        out.append(cls)
    return out


def reduce_label(label, has_inverse):
    bar = label.startswith("~")
    base = label[1:] if bar else label
    base = re.sub(r"\(\d+\)$", "", base)
    m = re.match(r"^([CS])(\d+)\^(-?\d+)$", base)
    if m:
        letter, order, j = m.group(1), int(m.group(2)), int(m.group(3))
        if has_inverse:
            j = abs(j)
        base = f"{letter}{order}" if j == 1 else f"{letter}{order}^{j}"
    return ("~" if bar else "") + base


def class_label(size, rep_label, has_inverse):
    return (str(size) if size > 1 else "") + reduce_label(rep_label, has_inverse)


# --------------------------------------------------------------------------
# Group context: everything the documents draw on
# --------------------------------------------------------------------------

FAMILY_NAMES = {
    "C": "Cn",
    "Cv": "Cnv",
    "Ch": "Cnh",
    "S": "S2n",
    "D": "Dn",
    "Dh": "Dnh",
    "Dd": "Dnd",
}
CONDUCTOR = {"T": 12, "O": 8, "I": 20}


def canonical_conductor(R):
    if R.name in CONDUCTOR:
        return CONDUCTOR[R.name]
    N = int(R.name[1:])
    if R.name.startswith("D"):
        return math.lcm(2 * N, 4)
    return N if N % 2 else 2 * N


class Context:
    def __init__(self, input_symbol):
        self.input_symbol = input_symbol
        self.G = ag.point_group(input_symbol)
        self.symbol = ag.ALIASES.get(input_symbol, input_symbol)
        self.family, self.n = ag.parse(self.symbol)
        self.R = self.G.rotation_group
        self.polyhedral = self.R.name in ("T", "O", "I")
        self.D = dg.double_group(self.G)
        self.N = self.G.order
        self.labels = ag.canonical_labels(self.G)
        self.labels_list = [self.labels[g] for g in self.G.elements]
        self.double_labels = self.labels_list + ["~" + s for s in self.labels_list]
        self.table = self.G.cayley()
        self.double_table = self.D.cayley()
        omega = self.D.factor_system()
        self.omega = [[omega[(a, b)] for b in self.G.elements] for a in self.G.elements]
        self.field = Field(canonical_conductor(self.R))
        self.classes = self._classes_of_G()
        self.double_classes = self._classes_of_double()
        self.irreps = build_irreps(self)

    # ---- classes
    def _classes_of_G(self):
        raw = classes_from_table(self.table)
        inv = inverses(self.table)
        cls_of = {g: c for c, members in enumerate(raw) for g in members}
        double_raw = classes_from_table(self.double_table)
        dcls_of = {q: c for c, members in enumerate(double_raw) for q in members}
        out = []
        for c, members in enumerate(raw):
            rep = members[0]
            has_inv = inv[rep] in members
            out.append(
                {
                    "index": c,
                    "label": class_label(len(members), self.labels_list[rep], has_inv),
                    "size": len(members),
                    "representative": rep,
                    "members": members,
                    "inverse": cls_of[inv[rep]],
                    "regular": dcls_of[rep] != dcls_of[rep + self.N],
                }
            )
        return out

    def _classes_of_double(self):
        raw = classes_from_table(self.double_table)
        inv = inverses(self.double_table)
        cls_of = {q: c for c, members in enumerate(raw) for q in members}
        ordered = []  # (base class index, representative, members, split, barred)
        for base in self.classes:
            g = base["representative"]
            lifted, barred = g, g + self.N
            split = cls_of[lifted] != cls_of[barred]
            for rep, is_barred in ((lifted, False), (barred, True)):
                if is_barred and not split:
                    break
                members = [rep] + [q for q in raw[cls_of[rep]] if q != rep]
                ordered.append((base["index"], rep, members, split, is_barred))
        index_of_rep = {rep: k for k, (_, rep, _, _, _) in enumerate(ordered)}
        out = []
        for k, (base, rep, members, split, is_barred) in enumerate(ordered):
            has_inv = inv[rep] in members
            if split:
                label = class_label(len(members), self.double_labels[rep], has_inv)
            else:
                half = len(members) // 2
                label = (
                    class_label(half, self.double_labels[rep], has_inv)
                    + "+"
                    + class_label(half, self.double_labels[rep + self.N], has_inv)
                )
            inv_rep = inv[rep]
            inv_class = next(
                kk for kk, (_, r2, m2, _, _) in enumerate(ordered) if inv_rep in m2
            )
            rec = {
                "index": k,
                "label": label,
                "size": len(members),
                "representative": rep,
                "members": members,
                "inverse": inv_class,
                "base_class": base,
                "split": split,
            }
            if split:
                rec["barred"] = is_barred
            out.append(rec)
        del index_of_rep
        return out

    # ---- records
    def element_records(self, double):
        classes = self.double_classes if double else self.classes
        cls_of = {g: c["index"] for c in classes for g in c["members"]}
        table = self.double_table if double else self.table
        inv = inverses(table)
        labels = self.double_labels if double else self.labels_list
        out = []
        for idx, label in enumerate(labels):
            base = idx % self.N
            rec = element_record(idx, label, self.G.elements[base][1], self.polyhedral)
            rec["inverse"] = inv[idx]
            rec["class"] = cls_of[idx]
            if double:
                rec["barred"] = idx >= self.N
                rec["base_index"] = base
            out.append(rec)
        return out

    def group_metadata(self):
        meta = {"symbol": self.symbol, "input_symbol": self.input_symbol}
        meta["family"] = self.symbol if self.polyhedral else FAMILY_NAMES[self.family]
        if not self.polyhedral:
            meta["n"] = self.n
        meta.update(order=self.N, type=self.G.kind, rotation_group=self.R.name)
        return meta


# --------------------------------------------------------------------------
# Irreps (docs/spec/irrep-labels.md, docs/spec/basis-conventions.md)
# --------------------------------------------------------------------------


class Irrep:
    def __init__(
        self,
        symbol,
        projective,
        chars,
        matrices=None,
        letter=None,
        subscript=None,
        suffix=None,
        component=None,
        partner_symbol=None,
    ):
        self.symbol = symbol
        self.projective = projective
        self.chars = (
            chars  # character on every element of G (lifted through the section)
        )
        self.matrices = matrices  # per element of G, or None (polyhedral prototype)
        self.letter = letter
        self.subscript = subscript
        self.suffix = suffix
        self.component = component
        self.partner_symbol = partner_symbol
        self.index = None
        self.partner = None
        self.frobenius_schur = None

    @property
    def dimension(self):
        return int(self.chars[0][0])

    def record(self):
        rec = {
            "index": self.index,
            "symbol": self.symbol,
            "dimension": self.dimension,
            "projective": self.projective,
            "letter": self.letter,
        }
        if self.component is not None:
            rec["component"] = self.component
            rec["partner"] = self.partner
        if self.subscript is not None:
            rec["subscript"] = self.subscript
        if self.suffix is not None:
            rec["suffix"] = self.suffix
        rec["frobenius_schur"] = self.frobenius_schur
        return rec


def _symbol_parts(symbol):
    m = re.match(r"^([12])?([ABETGHI])(\d+|\d*[13579]/2)?(g|u|'|'')?$", symbol)
    assert m, symbol
    comp, letter, sub, suf = m.groups()
    if sub is not None:
        sub = sub if "/" in sub else int(sub)
    return (int(comp) if comp else None), letter, sub, suf


def make_irrep(ctx, symbol, projective, chars, matrices=None, partner_symbol=None):
    comp, letter, sub, suf = _symbol_parts(symbol)
    return Irrep(
        symbol, projective, chars, matrices, letter, sub, suf, comp, partner_symbol
    )


LETTER_RANK = {L: k for k, L in enumerate("ABETGHI")}
SUFFIX_RANK = {None: 0, "g": 0, "u": 1, "'": 0, "''": 1}


def order_key(ir):
    sub = ir.subscript
    if sub is None:
        subv = -1.0
    elif isinstance(sub, str):
        subv = int(sub.split("/")[0]) / 2
    else:
        subv = float(sub)
    return (
        ir.projective,
        SUFFIX_RANK[ir.suffix],
        LETTER_RANK[ir.letter],
        subv,
        ir.component or 0,
    )


def with_parity(ctx, on_R, suffix_sign):
    """Lift a function on R (or on lifted R-dagger) to G: D((r, p)) = sign^(p<0) D(r)."""
    return [
        on_R(r) if (p == 1 or suffix_sign == 1) else neg_matrix(ctx.field, on_R(r))
        for (r, p) in ctx.G.elements
    ]


def neg_matrix(F, A):
    return [[F.neg(a) for a in row] for row in A]


def dihedral_irreps(ctx):
    """All irreps of G with R = D_N, from the canonical-basis formulas."""
    F = ctx.field
    N = int(ctx.R.name[1:])
    M = F.conductor
    w, u, i1 = F.zeta(M // N), F.zeta(M // (2 * N)), F.i()
    one = mat_identity(F, 1)

    def rep(Da, Db, lifted):  # r = (k, s) in D_N, or (k~, s) in Dic_N when lifted
        def D(r):
            k, s = r
            return mat_mul(F, mat_power(F, Da, k), mat_power(F, Db, s))

        if lifted:
            return lambda r: D(ctx.D.cover.section(r))
        return D

    reps = []  # (symbol, projective, function on R)
    reps.append(("A1", False, rep(one, one, False)))
    reps.append(("A2", False, rep(one, neg_matrix(F, one), False)))
    if N % 2 == 0:
        reps.append(("B1", False, rep(neg_matrix(F, one), one, False)))
        reps.append(("B2", False, rep(neg_matrix(F, one), neg_matrix(F, one), False)))
    for k in range(1, (N + 1) // 2):
        if 2 * k == N:
            break
        Da = [[F.power(w, k), F.zero()], [F.zero(), F.power(w, (-k) % N)]]
        Db = [[F.zero(), F.rational(1)], [F.rational(1), F.zero()]]
        reps.append(("E" if N <= 4 else f"E{k}", False, rep(Da, Db, False)))
    for k in range(1, N, 2):
        Da = [[F.power(u, k), F.zero()], [F.zero(), F.power(u, (-k) % (2 * N))]]
        Db = [[F.zero(), i1], [i1, F.zero()]]
        reps.append((f"E{k}/2", True, rep(Da, Db, True)))
    if N % 2 == 1:
        reps.append((f"1E{N}/2", True, rep(neg_matrix(F, one), [[i1]], True)))
        reps.append((f"2E{N}/2", True, rep(neg_matrix(F, one), [[F.neg(i1)]], True)))

    suffixes = [("", 1)] if ctx.G.kind != "centro" else [("g", 1), ("u", -1)]
    out = []
    for symbol, projective, D in reps:
        for suffix, sign in suffixes:
            mats = with_parity(ctx, D, sign)
            chars = [mat_trace(F, A) for A in mats]
            partner = None
            if symbol.startswith(("1", "2")):
                partner = ("2" if symbol[0] == "1" else "1") + symbol[1:] + suffix
            out.append(
                make_irrep(ctx, symbol + suffix, projective, chars, mats, partner)
            )
    return out


def polyhedral_irreps(ctx):
    """Characters of every irrep of O(h) and I(h) from the quaternion w-coordinate."""
    F = ctx.field
    kind = ctx.R.name
    alpha = F.sqrt2() if kind == "O" else F.phi()

    def w_of(r, conj=False):
        a, b = ctx.D.cover.section(r)[0]  # w = a + b*alpha
        if conj:  # Galois conjugate phi -> 1 - phi
            a, b = a + b, -b
        return F.add(F.rational(a), F.scale(alpha, b))

    def U(w, deg):  # Chebyshev U_deg(w): character of D^(deg/2) at cos(theta/2) = w
        coeffs = {
            1: [(2, 1)],
            2: [(4, 2), (-1, 0)],
            3: [(8, 3), (-4, 1)],
            4: [(16, 4), (-12, 2), (1, 0)],
            5: [(32, 5), (-32, 3), (6, 1)],
            6: [(64, 6), (-80, 4), (24, 2), (-1, 0)],
        }[deg]
        return _sum(F, (F.scale(F.power(w, e), c) for c, e in coeffs))

    def perm_sign(p):
        seen, sign = set(), 1
        for s in range(len(p)):
            if s in seen:
                continue
            length, x = 0, s
            while x not in seen:
                seen.add(x)
                x = p[x]
                length += 1
            sign *= (-1) ** (length - 1)
        return sign

    reps = []  # (symbol, projective, function r -> character)
    if kind == "O":
        A2 = lambda r: F.rational(perm_sign(r))
        T1 = lambda r: U(w_of(r), 2)
        T2 = lambda r: F.mul(A2(r), T1(r))
        E = lambda r: F.sub(
            F.sub(F.sub(F.mul(T1(r), T1(r)), F.rational(1)), T1(r)), T2(r)
        )
        reps += [
            ("A1", False, lambda r: F.rational(1)),
            ("A2", False, A2),
            ("E", False, E),
            ("T1", False, T1),
            ("T2", False, T2),
        ]
        E12 = lambda r: U(w_of(r), 1)
        reps += [
            ("E1/2", True, E12),
            ("E5/2", True, lambda r: F.mul(A2(r), E12(r))),
            ("G3/2", True, lambda r: U(w_of(r), 3)),
        ]
    else:
        T1 = lambda r: U(w_of(r), 2)
        T2 = lambda r: U(w_of(r, True), 2)
        reps += [
            ("A", False, lambda r: F.rational(1)),
            ("T1", False, T1),
            ("T2", False, T2),
            ("G", False, lambda r: F.sub(U(w_of(r), 6), T2(r))),
            ("H", False, lambda r: U(w_of(r), 4)),
        ]
        reps += [
            ("E1/2", True, lambda r: U(w_of(r), 1)),
            ("E7/2", True, lambda r: U(w_of(r, True), 1)),
            ("G3/2", True, lambda r: U(w_of(r), 3)),
            ("I5/2", True, lambda r: U(w_of(r), 5)),
        ]

    suffixes = [("", 1)] if ctx.G.kind != "centro" else [("g", 1), ("u", -1)]
    out = []
    for symbol, projective, chi in reps:
        for suffix, sign in suffixes:
            chars = [
                chi(r) if (p == 1 or sign == 1) else F.neg(chi(r))
                for (r, p) in ctx.G.elements
            ]
            out.append(make_irrep(ctx, symbol + suffix, projective, chars))
    return out


def build_irreps(ctx):
    if ctx.polyhedral:
        assert ctx.R.name in ("O", "I"), (
            "prototype covers O, Oh, I, Ih among the polyhedral groups"
        )
        irreps = polyhedral_irreps(ctx)
    else:
        assert ctx.R.name.startswith("D"), (
            "prototype covers the dihedral-based axial families"
        )
        irreps = dihedral_irreps(ctx)
    irreps.sort(key=order_key)
    by_symbol = {}
    for k, ir in enumerate(irreps):
        ir.index = k
        by_symbol[ir.symbol] = k
    F = ctx.field
    for ir in irreps:
        if ir.partner_symbol:
            ir.partner = by_symbol[ir.partner_symbol]
        # Frobenius-Schur: (1/|G|) sum_g omega(g, g) chi(g^2)
        total = F.zero()
        for a, g in enumerate(ctx.G.elements):
            sq = ctx.table[a][a]
            term = ir.chars[sq]
            if ir.projective and ctx.omega[a][a] == -1:
                term = F.neg(term)
            total = F.add(total, term)
        total = F.scale(total, Fr(1, ctx.N))
        assert all(c == 0 for c in total[1:]) and total[0].denominator == 1, (
            ir.symbol,
            total,
        )
        ir.frobenius_schur = int(total[0])
    return irreps


def character_on(ctx, ir, double_index):
    """Character of an irrep at an element of G (index < N) or G-dagger (any index)."""
    g = double_index % ctx.N
    chi = ir.chars[g]
    if double_index >= ctx.N and ir.projective:
        chi = ctx.field.neg(chi)
    return chi


def select_irreps(ctx, which):
    return [
        ir
        for ir in ctx.irreps
        if which == "all" or (which == "projective") == ir.projective
    ]


# --------------------------------------------------------------------------
# Documents (the layout under discussion)
# --------------------------------------------------------------------------

TABLE_ORDER = ["elements", "cayley", "factor", "characters", "matrices"]


def value_forms(values):
    forms = []
    for v in values.split(","):
        if v == "both":
            forms += ["exact", "numeric"]
        else:
            forms.append(v)
    return [f for f in ["exact", "numeric", "pretty"] if f in forms]


def elements_table(ctx, double):
    return {
        "group": "G_double" if double else "G",
        "items": ctx.element_records(double),
        "classes": ctx.double_classes if double else ctx.classes,
    }


def cayley_table(ctx, double):
    return {
        "group": "G_double" if double else "G",
        "table": ctx.double_table if double else ctx.table,
    }


def factor_table(ctx):
    return {"group": "G", "table": ctx.omega}


def characters_table(ctx, which, forms):
    F = ctx.field
    irreps = select_irreps(ctx, which)
    if which == "ordinary":
        class_group, classes = "G", ctx.classes
    else:
        class_group, classes = "G_double", ctx.double_classes
    rows = [
        [character_on(ctx, ir, c["representative"]) for c in classes] for ir in irreps
    ]
    table = {
        "class_group": class_group,
        "irreps": [ir.record() for ir in irreps],
        "classes": classes,
    }
    if "exact" in forms:
        table["exact"] = [[F.exact(x) for x in row] for row in rows]
    if "numeric" in forms:
        table["numeric"] = [[F.numeric(x) for x in row] for row in rows]
    return table


def matrices_table(ctx, which, forms):
    F = ctx.field
    entries = []
    for ir in select_irreps(ctx, which):
        rec = ir.record()
        if ir.matrices is None:
            raise SystemExit(
                f"prototype: no irrep matrices for {ctx.symbol} yet (ticket #16)"
            )
        if "exact" in forms:
            rec["exact"] = [mat_exact(F, A) for A in ir.matrices]
        if "numeric" in forms:
            rec["numeric"] = [mat_numeric(F, A) for A in ir.matrices]
        entries.append(rec)
    return {"irreps": entries}


def document(ctx, tables, double=False, irreps="all", values="both"):
    forms = value_forms(values)
    tables = [t for t in TABLE_ORDER if t in tables]
    if "elements" not in tables:  # every document embeds the element list (proposal)
        tables = ["elements"] + tables
    doc = {
        "schema_version": 1,
        "generator": {"name": "pointgroup-tables", "version": PROTOTYPE_VERSION},
        "group": ctx.group_metadata(),
        "tables": tables,
        "double": double,
    }
    if "characters" in tables or "matrices" in tables:
        doc["irreps"] = irreps
    doc["values"] = forms
    doc["field"] = ctx.field.metadata()
    doc["basis_convention"] = {"name": "canonical"}
    for t in tables:
        if t == "elements":
            doc[t] = elements_table(ctx, double)
        elif t == "cayley":
            doc[t] = cayley_table(ctx, double)
        elif t == "factor":
            doc[t] = factor_table(ctx)
        elif t == "characters":
            doc[t] = characters_table(ctx, irreps, forms)
        elif t == "matrices":
            doc[t] = matrices_table(ctx, irreps, forms)
    return doc


# --------------------------------------------------------------------------
# JSON printing: objects one key per line, leaf arrays inline
# --------------------------------------------------------------------------


def _inline(x):
    return json.dumps(x, separators=(", ", ": "))


def _is_leafy(x):
    if isinstance(x, dict):
        return False
    if isinstance(x, list):
        return all(_is_leafy(y) for y in x)
    return True


def dumps(x, indent=0, width=100):
    pad = " " * indent
    if isinstance(x, dict):
        if not x:
            return "{}"
        items = [f'{pad}  "{k}": {dumps(v, indent + 2, width)}' for k, v in x.items()]
        return "{\n" + ",\n".join(items) + f"\n{pad}}}"
    if isinstance(x, list):
        if not x:
            return "[]"
        if _is_leafy(x):
            s = _inline(x)
            if len(s) + indent <= width:
                return s
        items = [f"{pad}  {dumps(v, indent + 2, width)}" for v in x]
        return "[\n" + ",\n".join(items) + f"\n{pad}]"
    return _inline(x)


# --------------------------------------------------------------------------
# Checks
# --------------------------------------------------------------------------


def selfcheck(symbols=("C3v", "D6h", "Oh", "Ih")):
    for symbol in symbols:
        ctx = Context(symbol)
        F = ctx.field
        irreps = ctx.irreps
        for proj in (False, True):
            block = [ir for ir in irreps if ir.projective == proj]
            assert sum(ir.dimension**2 for ir in block) == ctx.N, (symbol, proj)
        # orthogonality over G-dagger (all irreps are irreps of G-dagger)
        classes = ctx.double_classes
        for a in irreps:
            for b in irreps:
                total = F.zero()
                for c in classes:
                    x = character_on(ctx, a, c["representative"])
                    y = F.conj(character_on(ctx, b, c["representative"]))
                    total = F.add(total, F.scale(F.mul(x, y), c["size"]))
                expected = F.rational(2 * ctx.N if a is b else 0)
                assert total == expected, (symbol, a.symbol, b.symbol, total)
        # D(g) D(h) = omega(g, h) D(gh) for the axial matrices
        for ir in irreps:
            if ir.matrices is None:
                continue
            for a in range(ctx.N):
                for b in range(ctx.N):
                    prod = mat_mul(F, ir.matrices[a], ir.matrices[b])
                    target = ir.matrices[ctx.table[a][b]]
                    if ir.projective and ctx.omega[a][b] == -1:
                        target = neg_matrix(F, target)
                    assert prod == target, (symbol, ir.symbol, a, b)
        nreg = sum(1 for c in ctx.classes if c["regular"])
        nproj = sum(1 for ir in irreps if ir.projective)
        assert nreg == nproj, (symbol, nreg, nproj)
        print(
            f"{symbol:4s} ok: |G|={ctx.N:3d} classes G/G†={len(ctx.classes)}/{len(classes)} "
            f"irreps={' '.join(ir.symbol for ir in irreps)}"
        )


# --------------------------------------------------------------------------
# Command line (mirrors docs/spec/cli.md for the table subcommands)
# --------------------------------------------------------------------------

SAMPLES = [  # (subcommand, group, options)
    ("elements", "C3v", {}),
    ("elements", "C3v", {"double": True}),
    ("cayley", "C3v", {}),
    ("factor", "C3v", {}),
    ("characters", "C3v", {}),
    ("characters", "C3v", {"irreps": "ordinary"}),
    ("matrices", "C3v", {}),
    ("all", "C3v", {}),
    ("elements", "D6h", {"double": True}),
    ("characters", "D6h", {}),
    ("matrices", "D6h", {"values": "exact"}),
    ("all", "D6h", {}),
    ("elements", "Oh", {"double": True}),
    ("characters", "Oh", {}),
    ("all", "Oh", {"tables": ["elements", "cayley", "factor", "characters"]}),
    ("elements", "Ih", {}),
    ("characters", "Ih", {}),
    ("all", "Ih", {"tables": ["elements", "cayley", "factor", "characters"]}),
]


def write_samples(out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)
    contexts = {}
    rows = []
    for sub, group, opts in SAMPLES:
        ctx = contexts.setdefault(group, Context(group))
        tables = opts.get("tables", TABLE_ORDER if sub == "all" else [sub])
        doc = document(
            ctx,
            tables,
            opts.get("double", False),
            opts.get("irreps", "all"),
            opts.get("values", "both"),
        )
        name = f"{group}_{sub}"
        if opts.get("double"):
            name += "_double"
        if opts.get("irreps"):
            name += "_" + str(opts["irreps"])
        if opts.get("values"):
            name += "_" + str(opts["values"])
        path = out_dir / f"{name}.json"
        text = dumps(doc) + "\n"
        path.write_text(text)
        rows.append((path.name, len(text)))
    for name, size in rows:
        print(f"{size:9d}  {name}")


def main(argv):
    p = argparse.ArgumentParser(prog="output_schema_prototype")
    p.add_argument("subcommand", choices=TABLE_ORDER + ["all", "selfcheck", "samples"])
    p.add_argument("group", nargs="?")
    p.add_argument("--double-group", action="store_true")
    p.add_argument("--irreps", choices=["ordinary", "projective", "all"], default="all")
    p.add_argument("--values", default="both")
    p.add_argument("-o", "--output")
    args = p.parse_args(argv)
    if args.subcommand == "selfcheck":
        return selfcheck()
    if args.subcommand == "samples":
        return write_samples(Path(__file__).resolve().parent / "samples")
    if not args.group:
        p.error("a group is required")
    ctx = Context(args.group)
    tables = TABLE_ORDER if args.subcommand == "all" else [args.subcommand]
    doc = document(ctx, tables, args.double_group, args.irreps, args.values)
    text = dumps(doc) + "\n"
    if args.output and args.output != "-":
        Path(args.output).write_text(text)
    else:
        sys.stdout.write(text)


if __name__ == "__main__":
    main(sys.argv[1:])
