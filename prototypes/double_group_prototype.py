"""Prototype: spin double groups and the spin factor system (wayfinder ticket #4).

Builds the double group G-dagger of every point group produced by
``abstract_group_prototype`` as the same subgroup construction applied to
R-dagger x {+1, -1}, where R-dagger is the binary cover of the rotation group:

* Z_2n for C_n (residues mod 2n),
* Dic_n for D_n (pairs (k mod 2n, s) with b^2 = e-bar),
* 2T, 2O, 2I as unit quaternions with coordinates in (1/2) Z[alpha],
  alpha = 1, sqrt2, phi, generated from fixed quaternions whose images are the
  permutation generators of T, O, I.

The canonical section s(r) is the *positive* preimage: residue in
(-n/2, n/2] for the axial covers, first non-zero quaternion coordinate
positive for the polyhedral covers. The factor system is
s(g) s(h) = omega(g, h) s(gh).

Run without arguments for the self-check over 43 groups (group axioms of
G-dagger, normalised cocycle, class counts), or pass symbols to print the
double-group element list and the factor table:

    uv run python prototypes/double_group_prototype.py
    uv run python prototypes/double_group_prototype.py C3v D3h T

Everything is integer arithmetic; see docs/spec/double-groups.md.
"""

from __future__ import annotations

import sys
from collections.abc import Callable
from dataclasses import dataclass
from fractions import Fraction as Fr
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import abstract_group_prototype as ag

# --------------------------------------------------------------------------
# The quadratic ring Z[alpha] with halves: value = a + b*alpha, alpha^2 = c*alpha + d
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Ring:
    name: str
    c: int
    d: int

    def mul(self, x, y):
        a, b = x
        e, f = y
        return (a * e + b * f * self.d, a * f + b * e + b * f * self.c)

    def sign(self, x):
        """Exact sign of a + b*alpha as a real number (alpha = 1, sqrt2 or phi)."""
        a, b = x
        if self.name == "one":
            v = a + b
            return (v > 0) - (v < 0)
        if self.name == "phi":  # a + b phi = ((2a + b) + b sqrt5) / 2
            a, d = 2 * a + b, 5
        else:
            d = self.d
        if a >= 0 and b >= 0:
            return 1 if (a, b) != (0, 0) else 0
        if a <= 0 and b <= 0:
            return -1
        return (1 if a * a > d * b * b else -1) * (1 if a > 0 else -1)


ONE = Ring("one", 0, 1)
SQRT2 = Ring("sqrt2", 0, 2)
PHI = Ring("phi", 1, 1)


def const(a, b=0):
    return (Fr(a), Fr(b))


ZERO = const(0)
QUAT_ONE = (const(1), ZERO, ZERO, ZERO)
QUAT_EBAR = (const(-1), ZERO, ZERO, ZERO)


def qmul(ring: Ring, p, q):
    w1, x1, y1, z1 = p
    w2, x2, y2, z2 = q
    m = ring.mul

    def add(*terms):
        return (sum(t[0] for t in terms), sum(t[1] for t in terms))

    def neg(t):
        return (-t[0], -t[1])

    return (
        add(m(w1, w2), neg(m(x1, x2)), neg(m(y1, y2)), neg(m(z1, z2))),
        add(m(w1, x2), m(x1, w2), m(y1, z2), neg(m(z1, y2))),
        add(m(w1, y2), neg(m(x1, z2)), m(y1, w2), m(z1, x2)),
        add(m(w1, z2), m(x1, y2), neg(m(y1, x2)), m(z1, w2)),
    )


def qsign(ring: Ring, q):
    for coord in q:
        s = ring.sign(coord)
        if s:
            return s
    raise ValueError("zero quaternion")


# --------------------------------------------------------------------------
# Binary rotation groups R-dagger with projection and canonical section
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class BinaryRotationGroup:
    name: str
    elements: tuple
    mul: Callable
    ebar: object
    proj: Callable  # R-dagger -> R
    section: Callable  # R -> R-dagger, the positive preimage


def cyclic_double(n: int) -> BinaryRotationGroup:
    N = 2 * n

    def section(k):
        return k % N if 2 * k <= n else (k - n) % N

    return BinaryRotationGroup(
        name=f"Z{N}",
        elements=tuple(range(N)),
        mul=lambda a, b: (a + b) % N,
        ebar=n,
        proj=lambda a: a % n,
        section=section,
    )


def dicyclic_double(n: int) -> BinaryRotationGroup:
    N = 2 * n

    def mul(a, b):
        k1, s1 = a
        k2, s2 = b
        return ((k1 + (k2 if s1 == 0 else -k2) + s1 * s2 * n) % N, (s1 + s2) % 2)

    def section(r):
        k, s = r
        return ((k if 2 * k <= n else k - n) % N, s)

    elements = tuple((k, s) for s in (0, 1) for k in range(N))
    return BinaryRotationGroup(
        name=f"Dic{n}",
        elements=elements,
        mul=mul,
        ebar=(n, 0),
        proj=lambda a: (a[0] % n, a[1]),
        section=section,
    )


# Fixed quaternion generators; their images are the permutation generators of
# abstract_group_prototype.polyhedral (docs/spec/double-groups.md).
half = Fr(1, 2)
GENERATORS = {
    "T": (
        ONE,
        [
            ((const(half), const(half), const(half), const(half)), (1, 2, 0, 3)),
            ((ZERO, ZERO, ZERO, const(1)), (1, 0, 3, 2)),
        ],
    ),
    "O": (
        SQRT2,
        [
            ((const(0, half), ZERO, ZERO, const(0, half)), (1, 2, 3, 0)),
            ((ZERO, ZERO, const(0, half), const(0, half)), (1, 0, 2, 3)),
        ],
    ),
    "I": (
        PHI,
        [
            # (phi, 1/phi, 0, 1)/2 and (1, 1/phi, phi, 0)/2 with 1/phi = phi - 1
            ((const(0, half), const(-half, half), ZERO, const(half)), (1, 2, 3, 4, 0)),
            ((const(half), const(-half, half), const(0, half), ZERO), (1, 2, 0, 3, 4)),
        ],
    ),
}


def polyhedral_double(kind: str) -> BinaryRotationGroup:
    ring, gens = GENERATORS[kind]
    proj = {QUAT_ONE: tuple(range(len(gens[0][1])))}
    elements = [QUAT_ONE]
    frontier = [QUAT_ONE]
    while frontier:  # breadth-first closure with the homomorphism tracked
        nxt = []
        for g in frontier:
            for q, perm in gens:
                gq = qmul(ring, g, q)
                image = ag._perm_mul(proj[g], perm)
                if gq not in proj:
                    proj[gq] = image
                    elements.append(gq)
                    nxt.append(gq)
                elif proj[gq] != image:
                    raise AssertionError(
                        f"generator images of 2{kind} are not a homomorphism"
                    )
        frontier = nxt
    expected = {"T": 24, "O": 48, "I": 120}[kind]
    assert len(elements) == expected, (kind, len(elements))
    for q in elements:  # pi(q) = pi(-q), fibres of size 2
        assert proj[q] == proj[qmul(ring, QUAT_EBAR, q)]
    assert len(set(proj.values())) == expected // 2
    positive = {}
    for q in elements:
        if qsign(ring, q) > 0:
            assert proj[q] not in positive
            positive[proj[q]] = q
    return BinaryRotationGroup(
        name=f"2{kind}",
        elements=tuple(elements),
        mul=lambda a, b: qmul(ring, a, b),
        ebar=QUAT_EBAR,
        proj=lambda q: proj[q],
        section=lambda r: positive[r],
    )


_CACHE = {}


def binary_cover(R: ag.RotationGroup) -> BinaryRotationGroup:
    if R.name not in _CACHE:
        if R.name in ("T", "O", "I"):
            _CACHE[R.name] = polyhedral_double(R.name)
        elif R.name.startswith("D"):
            _CACHE[R.name] = dicyclic_double(int(R.name[1:]))
        else:
            _CACHE[R.name] = cyclic_double(int(R.name[1:]))
    return _CACHE[R.name]


# --------------------------------------------------------------------------
# Double group G-dagger and factor system
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class DoubleGroup:
    group: ag.PointGroup
    cover: BinaryRotationGroup
    elements: tuple  # lifted elements in canonical order, then barred elements

    @property
    def order(self):
        return len(self.elements)

    def mul(self, a, b):
        return (self.cover.mul(a[0], b[0]), a[1] * b[1])

    def section(self, g):
        r, p = g
        return (self.cover.section(r), p)

    def bar(self, q):
        return (self.cover.mul(self.cover.ebar, q[0]), q[1])

    def index(self, q):
        return self.elements.index(q)

    def cayley(self):
        idx = {q: i for i, q in enumerate(self.elements)}
        return [[idx[self.mul(a, b)] for b in self.elements] for a in self.elements]

    def factor_system(self):
        G = self.group
        omega = {}
        for a in G.elements:
            for b in G.elements:
                prod = self.mul(self.section(a), self.section(b))
                lift = self.section(G.mul(a, b))
                if prod == lift:
                    omega[(a, b)] = 1
                elif prod == self.bar(lift):
                    omega[(a, b)] = -1
                else:
                    raise AssertionError(
                        "product of lifts is not a preimage of the product"
                    )
        return omega

    def labels(self):
        base = ag.canonical_labels(self.group)
        n = self.group.order
        out = {}
        for i, q in enumerate(self.elements):
            g = self.group.elements[i % n]
            out[q] = base[g] if i < n else "~" + base[g]
        return out


def double_group(G: ag.PointGroup) -> DoubleGroup:
    cover = binary_cover(G.rotation_group)
    lifted = [(cover.section(r), p) for (r, p) in G.elements]
    barred = [(cover.mul(cover.ebar, q), p) for (q, p) in lifted]
    return DoubleGroup(G, cover, tuple(lifted + barred))


# --------------------------------------------------------------------------
# Checks
# --------------------------------------------------------------------------


def check_double_group(D: DoubleGroup):
    table = D.cayley()
    n = len(table)
    assert len(set(D.elements)) == n == 2 * D.group.order
    assert all(table[0][j] == j and table[j][0] == j for j in range(n))
    assert all(0 in table[i] for i in range(n))
    if n <= 96:
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    assert table[table[i][j]][k] == table[i][table[j][k]]
    return table


def check_factor_system(G: ag.PointGroup, omega):
    E = G.identity()
    assert set(omega.values()) <= {1, -1}
    for g in G.elements:
        assert omega[(E, g)] == 1 and omega[(g, E)] == 1
        assert omega[(g, G.inv(g))] == omega[(G.inv(g), g)]
    els = G.elements if G.order <= 48 else G.elements[::5]
    for g in els:
        for h in els:
            for k in els:
                assert (
                    omega[(g, h)] * omega[(G.mul(g, h), k)]
                    == omega[(g, G.mul(h, k))] * omega[(h, k)]
                )


def class_count(table):
    n = len(table)
    inv = [table[i].index(0) for i in range(n)]
    seen, classes = set(), 0
    for g in range(n):
        if g in seen:
            continue
        seen |= {table[table[h][g]][inv[h]] for h in range(n)}
        classes += 1
    return classes


# (classes of G, classes of G-dagger); Dic_n has n + 3 classes.
KNOWN = {
    "C1": (1, 2),
    "Ci": (2, 4),
    "C2": (2, 4),
    "Cs": (2, 4),
    "C2v": (4, 5),
    "C2h": (4, 8),
    "D2": (4, 5),
    "D2h": (8, 10),
    "C3": (3, 6),
    "C3v": (3, 6),
    "D3": (3, 6),
    "D3d": (6, 12),
    "C3h": (6, 12),
    "D3h": (6, 9),
    "S6": (6, 12),
    "C4": (4, 8),
    "C4v": (5, 7),
    "D4": (5, 7),
    "D2d": (5, 7),
    "C4h": (8, 16),
    "D4h": (10, 14),
    "S4": (4, 8),
    "C6": (6, 12),
    "C6v": (6, 9),
    "D6": (6, 9),
    "C6h": (12, 24),
    "D6h": (12, 18),
    "C5": (5, 10),
    "C5v": (4, 8),
    "D5h": (8, 13),
    "D5d": (8, 16),
    "D8d": (11, 19),
    "C7v": (5, 10),
    "D7h": (10, 17),
    "C8h": (16, 32),
    "T": (4, 7),
    "Td": (5, 8),
    "Th": (8, 14),
    "O": (5, 8),
    "Oh": (10, 16),
    "I": (5, 9),
    "Ih": (10, 18),
}


def selfcheck():
    failures = 0
    for symbol, (nc, ncd) in KNOWN.items():
        G = ag.point_group(symbol)
        D = double_group(G)
        table = check_double_group(D)
        omega = D.factor_system()
        check_factor_system(G, omega)
        got = (len(ag.conjugacy_classes(G)), class_count(table))
        ok = got == (nc, ncd)
        failures += not ok
        neg = sum(v == -1 for v in omega.values())
        print(
            f"{symbol:4s} |G|={G.order:3d} {D.cover.name:6s} classes {got[0]:2d}/{got[1]:2d} "
            f"projective irreps {got[1] - got[0]:2d} omega=-1 on {neg:4d}/{G.order**2:5d} "
            f"{'ok' if ok else 'MISMATCH expected ' + str((nc, ncd))}"
        )
    print(f"{len(KNOWN) - failures}/{len(KNOWN)} groups match")
    return failures == 0


def show(symbol: str):
    G = ag.point_group(symbol)
    D = double_group(G)
    labels = D.labels()
    print(f"{G.symbol}-dagger: order {D.order}, cover {D.cover.name}")
    for i, q in enumerate(D.elements):
        print(f"  {i:3d} {labels[q]:14s} {q}")
    omega = D.factor_system()
    base = ag.canonical_labels(G)
    width = max(len(v) for v in base.values())
    print("factor system omega(row, column):")
    print(" " * (width + 1) + " ".join(f"{base[h]:>{width}s}" for h in G.elements))
    for g in G.elements:
        row = " ".join(f"{omega[(g, h)]:>{width}d}" for h in G.elements)
        print(f"{base[g]:{width}s} {row}")


if __name__ == "__main__":
    if len(sys.argv) == 1:
        sys.exit(0 if selfcheck() else 1)
    for sym in sys.argv[1:]:
        show(sym)
        print()
