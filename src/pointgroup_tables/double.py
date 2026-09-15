"""The spin double group ``G-dagger`` and the spin factor system
(docs/spec/double-groups.md, ADR-0003).

``G-dagger`` is the subgroup construction of the abstract group applied to
``R-dagger x {+1, -1}``, where ``R-dagger`` is the binary cover of the rotation
group: ``Z_2n`` for ``C_n``, ``Dic_n`` for ``D_n``, and the binary polyhedral
groups ``2T``, ``2O``, ``2I`` as unit quaternions with coordinates in
``(1/2) Z[alpha]``, ``alpha = 1, sqrt2, phi``. Everything is integer
arithmetic. The canonical section lifts every rotation to its *positive*
preimage (residue in ``(-n/2, n/2]``, or first non-zero quaternion coordinate
positive); the inversion lifts trivially.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from functools import cache
from typing import Any

from . import rotation
from .abstract import Element, PointGroup

# ---- the quadratic ring Z[alpha] with alpha^2 = c alpha + d ------------------

Coordinate = tuple[int, int]  # (a, b) meaning (a + b alpha) / 2
Quaternion = tuple[Coordinate, Coordinate, Coordinate, Coordinate]


@dataclass(frozen=True)
class Ring:
    """``Z[alpha]`` with ``alpha^2 = c alpha + d`` and an exact sign of ``a + b alpha``."""

    name: str
    c: int
    d: int

    def sign(self, a: int, b: int) -> int:
        """The sign of ``a + b alpha`` as a real number (``alpha = 1, sqrt2 or phi``)."""
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

ZERO: Coordinate = (0, 0)
QUATERNION_ONE: Quaternion = ((2, 0), ZERO, ZERO, ZERO)
QUATERNION_EBAR: Quaternion = ((-2, 0), ZERO, ZERO, ZERO)


def quaternion_mul(ring: Ring, p: Quaternion, q: Quaternion) -> Quaternion:
    """The Hamilton product; every coordinate is ``(a + b alpha) / 2``."""
    c, d = ring.c, ring.d

    def product(x: Coordinate, y: Coordinate) -> tuple[int, int]:
        # (a + b alpha)(e + f alpha) = (ae + bf d) + (af + be + bf c) alpha
        a, b = x
        e, f = y
        return (a * e + b * f * d, a * f + b * e + b * f * c)

    def combine(*terms: tuple[int, tuple[int, int]]) -> Coordinate:
        # the sum of signed products carries the denominator 4; the result has
        # denominator 2 by closure of the group, so both parts are even
        a = sum(sign * term[0] for sign, term in terms)
        b = sum(sign * term[1] for sign, term in terms)
        if a % 2 or b % 2:
            raise ArithmeticError("the product left the half-integer lattice")
        return (a // 2, b // 2)

    w1, x1, y1, z1 = p
    w2, x2, y2, z2 = q
    return (
        combine((1, product(w1, w2)), (-1, product(x1, x2)), (-1, product(y1, y2)), (-1, product(z1, z2))),
        combine((1, product(w1, x2)), (1, product(x1, w2)), (1, product(y1, z2)), (-1, product(z1, y2))),
        combine((1, product(w1, y2)), (-1, product(x1, z2)), (1, product(y1, w2)), (1, product(z1, x2))),
        combine((1, product(w1, z2)), (1, product(x1, y2)), (-1, product(y1, x2)), (1, product(z1, w2))),
    )  # fmt: skip


def quaternion_sign(ring: Ring, q: Quaternion) -> int:
    """The sign of the first non-zero coordinate in the order ``(w, x, y, z)``."""
    for a, b in q:
        sign = ring.sign(a, b)
        if sign:
            return sign
    raise ValueError("the zero quaternion has no sign")


# ---- the binary rotation groups -----------------------------------------------


@dataclass(frozen=True)
class BinaryRotationGroup:
    """``R-dagger`` with its projection onto ``R`` and the canonical section."""

    name: str
    elements: tuple[Any, ...]
    mul: Callable[[Any, Any], Any]
    ebar: Any
    project: Callable[[Any], Any]  # R-dagger -> R
    section: Callable[[Any], Any]  # R -> R-dagger, the positive preimage

    @property
    def identity(self) -> Any:
        return self.elements[0]


def cyclic_cover(n: int) -> BinaryRotationGroup:
    """``Z_2n`` over ``C_n``: residues mod ``2n``, ``ebar = n``."""
    N = 2 * n

    def section(k: int) -> int:
        return k % N if 2 * k <= n else (k - n) % N

    return BinaryRotationGroup(
        name=f"Z{N}",
        elements=tuple(range(N)),
        mul=lambda a, b: (a + b) % N,
        ebar=n,
        project=lambda a: a % n,
        section=section,
    )


def dicyclic_cover(n: int) -> BinaryRotationGroup:
    """``Dic_n`` over ``D_n``: pairs ``(k mod 2n, s)`` meaning ``a^k b^s`` with ``b^2 = ebar``."""
    N = 2 * n

    def mul(a: tuple[int, int], b: tuple[int, int]) -> tuple[int, int]:
        k1, s1 = a
        k2, s2 = b
        return ((k1 + (k2 if s1 == 0 else -k2) + s1 * s2 * n) % N, (s1 + s2) % 2)

    def section(r: tuple[int, int]) -> tuple[int, int]:
        k, s = r
        return ((k if 2 * k <= n else k - n) % N, s)

    elements = tuple((k, s) for s in (0, 1) for k in range(N))
    return BinaryRotationGroup(
        name=f"Dic{n}",
        elements=elements,
        mul=mul,
        ebar=(n, 0),
        project=lambda a: (a[0] % n, a[1]),
        section=section,
    )


# Fixed quaternion generators whose images are the permutation generators of
# rotation.polyhedral; part of the specification (docs/spec/double-groups.md).
_HALF: Coordinate = (1, 0)
_HALF_ALPHA: Coordinate = (0, 1)
_QUATERNION_GENERATORS: dict[
    str, tuple[Ring, int, list[tuple[Quaternion, rotation.Perm]]]
] = {  # ring, |R-dagger|, generators with their permutation images
    "T": (
        ONE,
        24,
        [
            ((_HALF, _HALF, _HALF, _HALF), (1, 2, 0, 3)),  # (1, 1, 1, 1) / 2
            ((ZERO, ZERO, ZERO, (2, 0)), (1, 0, 3, 2)),  # (0, 0, 0, 1)
        ],
    ),
    "O": (
        SQRT2,
        48,
        [
            ((_HALF_ALPHA, ZERO, ZERO, _HALF_ALPHA), (1, 2, 3, 0)),  # (sqrt2 / 2)(1, 0, 0, 1)
            ((ZERO, ZERO, _HALF_ALPHA, _HALF_ALPHA), (1, 0, 2, 3)),  # (sqrt2 / 2)(0, 0, 1, 1)
        ],
    ),
    "I": (
        PHI,
        120,
        [
            # (phi, 1/phi, 0, 1) / 2 and (1, 1/phi, phi, 0) / 2 with 1/phi = phi - 1
            ((_HALF_ALPHA, (-1, 1), ZERO, _HALF), (1, 2, 3, 4, 0)),
            ((_HALF, (-1, 1), _HALF_ALPHA, ZERO), (1, 2, 0, 3, 4)),
        ],
    ),
}  # fmt: skip


def polyhedral_cover(name: str) -> BinaryRotationGroup:
    """``2T``, ``2O`` or ``2I`` by breadth-first closure, with the projection tracked."""
    ring, expected, generators = _QUATERNION_GENERATORS[name]
    degree = len(generators[0][1])
    projection: dict[Quaternion, rotation.Perm] = {QUATERNION_ONE: tuple(range(degree))}
    elements = [QUATERNION_ONE]
    frontier = [QUATERNION_ONE]
    while frontier:
        following = []
        for g in frontier:
            for q, perm in generators:
                gq = quaternion_mul(ring, g, q)
                image = rotation.perm_mul(projection[g], perm)
                if gq not in projection:
                    projection[gq] = image
                    elements.append(gq)
                    following.append(gq)
                elif projection[gq] != image:
                    raise ArithmeticError(
                        f"the generator images of 2{name} are not a homomorphism"
                    )
        frontier = following
    if len(elements) != expected:
        raise ArithmeticError(
            f"2{name} closed with {len(elements)} elements, not {expected}"
        )
    positive: dict[rotation.Perm, Quaternion] = {}
    for q in elements:
        if projection[q] != projection[quaternion_mul(ring, QUATERNION_EBAR, q)]:
            raise ArithmeticError(
                f"the projection of 2{name} is not constant on fibres"
            )
        if quaternion_sign(ring, q) > 0:
            if projection[q] in positive:
                raise ArithmeticError(f"two positive preimages in 2{name}")
            positive[projection[q]] = q
    if len(positive) != expected // 2:
        raise ArithmeticError(f"the projection of 2{name} is not onto {name}")
    return BinaryRotationGroup(
        name=f"2{name}",
        elements=tuple(elements),
        mul=lambda a, b: quaternion_mul(ring, a, b),
        ebar=QUATERNION_EBAR,
        project=lambda q: projection[q],
        section=lambda r: positive[r],
    )


@cache
def _cover(name: str, axis_order: int | None) -> BinaryRotationGroup:
    if axis_order is None:
        return polyhedral_cover(name)
    if name.startswith("D"):
        return dicyclic_cover(axis_order)
    return cyclic_cover(axis_order)


def binary_cover(R: rotation.RotationGroup) -> BinaryRotationGroup:
    """The binary cover ``R-dagger`` of a rotation group, built once per ``R``."""
    return _cover(R.name, R.axis_order)


# ---- the double group and the factor system ------------------------------------

DoubleElement = tuple[Any, int]


@dataclass(frozen=True)
class DoubleGroup:
    """``G-dagger``: the lifted elements in canonical order, then the barred ones."""

    group: PointGroup
    cover: BinaryRotationGroup
    elements: tuple[DoubleElement, ...]

    @property
    def order(self) -> int:
        return len(self.elements)

    @property
    def identity(self) -> DoubleElement:
        return self.elements[0]

    def mul(self, a: DoubleElement, b: DoubleElement) -> DoubleElement:
        return (self.cover.mul(a[0], b[0]), a[1] * b[1])

    def inv(self, a: DoubleElement) -> DoubleElement:
        q, parity = a
        power, inverse = q, self.cover.identity
        while power != self.cover.identity:
            inverse = power
            power = self.cover.mul(power, q)
        return (inverse, parity)

    def section(self, g: Element) -> DoubleElement:
        return (self.cover.section(g[0]), g[1])

    def cayley(self) -> list[list[int]]:
        """``table[i][j]`` is the index of the product of elements ``i`` and ``j`` of ``G-dagger``."""
        index = {q: i for i, q in enumerate(self.elements)}
        return [[index[self.mul(a, b)] for b in self.elements] for a in self.elements]

    def factor_system(self) -> list[list[int]]:
        """``omega(g, h)`` in ``{1, -1}`` with ``s(g) s(h) = omega(g, h) s(gh)``, over ``G``."""
        G = self.group
        N = G.order
        index = {q: i for i, q in enumerate(self.elements)}
        lifted = self.elements[:N]
        table = []
        for a in lifted:
            row = []
            for b in lifted:
                product = index[self.mul(a, b)]
                row.append(1 if product < N else -1)
            table.append(row)
        return table


def double_group(G: PointGroup) -> DoubleGroup:
    """The double group of ``G`` with the canonical section."""
    cover = binary_cover(G.rotation_group)
    lifted = [(cover.section(r), p) for r, p in G.elements]
    barred = [(cover.mul(cover.ebar, q), p) for q, p in lifted]
    return DoubleGroup(G, cover, tuple(lifted + barred))
