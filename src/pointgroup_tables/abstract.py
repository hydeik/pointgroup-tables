"""Point groups as abstract groups ``G <= R x {+1, -1}`` (docs/spec/abstract-groups.md).

Every element is the pair ``(r, parity)``. The group is proper (``{(r, +1)}``),
centro (``R x {+1, -1}``) or twisted (``{(r, chi(r))}`` with ``chi = +1`` on an
index-2 subgroup ``H``). Products and inverses are componentwise. The element
order is the canonical one: proper elements first, then improper, each block
in the order of ``R``.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from typing import Any

from . import rotation
from .symbols import GroupSpec

Element = tuple[Any, int]

PROPER, CENTRO, TWISTED = "proper", "centro", "twisted"


@dataclass(frozen=True)
class PointGroup:
    """The abstract group of a point group, with its rotation group and type."""

    spec: GroupSpec
    rotation_group: rotation.RotationGroup
    type: str
    elements: tuple[Element, ...]

    @property
    def order(self) -> int:
        return len(self.elements)

    @property
    def identity(self) -> Element:
        return self.elements[0]

    def mul(self, a: Element, b: Element) -> Element:
        return (self.rotation_group.mul(a[0], b[0]), a[1] * b[1])

    def inv(self, a: Element) -> Element:
        return (self.rotation_group.inv(a[0]), a[1])

    def cayley(self) -> list[list[int]]:
        """``table[i][j]`` is the index of ``g_i g_j``."""
        index = {g: i for i, g in enumerate(self.elements)}
        return [[index[self.mul(a, b)] for b in self.elements] for a in self.elements]


def _proper(spec: GroupSpec, R: rotation.RotationGroup) -> PointGroup:
    return PointGroup(spec, R, PROPER, tuple((r, +1) for r in R.elements))


def _centro(spec: GroupSpec, R: rotation.RotationGroup) -> PointGroup:
    elements = tuple((r, +1) for r in R.elements) + tuple((r, -1) for r in R.elements)
    return PointGroup(spec, R, CENTRO, elements)


def _twisted(
    spec: GroupSpec, R: rotation.RotationGroup, in_h: Callable[[Any], bool]
) -> PointGroup:
    elements = tuple((r, +1) for r in R.elements if in_h(r)) + tuple(
        (r, -1) for r in R.elements if not in_h(r)
    )
    return PointGroup(spec, R, TWISTED, elements)


def _even_k(r: Any) -> bool:
    """The index-2 subgroup ``{k even}`` of ``C_2n`` or ``D_2n``."""
    k = r if isinstance(r, int) else r[0]
    return k % 2 == 0


def _no_flip(r: tuple[int, int]) -> bool:
    """The index-2 subgroup ``C_n = {(k, 0)}`` of ``D_n``."""
    return r[1] == 0


def _cnh(spec: GroupSpec, n: int) -> PointGroup:
    if n % 2 == 0:
        return _centro(spec, rotation.cyclic(n))
    return _twisted(spec, rotation.cyclic(2 * n), _even_k)


def _s2n(spec: GroupSpec, n: int) -> PointGroup:
    """``S2n`` with subscript ``n``: ``R = C_{n/2}`` for odd ``n/2``, else ``C_n`` twisted."""
    if (n // 2) % 2:
        return _centro(spec, rotation.cyclic(n // 2))
    return _twisted(spec, rotation.cyclic(n), _even_k)


def _dnh(spec: GroupSpec, n: int) -> PointGroup:
    if n % 2 == 0:
        return _centro(spec, rotation.dihedral(n))
    return _twisted(spec, rotation.dihedral(2 * n), _even_k)


def _dnd(spec: GroupSpec, n: int) -> PointGroup:
    if n % 2:
        return _centro(spec, rotation.dihedral(n))
    return _twisted(spec, rotation.dihedral(2 * n), _even_k)


# the family table of docs/spec/abstract-groups.md
AXIAL_FAMILY_TABLE: dict[str, Callable[[GroupSpec, int], PointGroup]] = {
    "Cn": lambda spec, n: _proper(spec, rotation.cyclic(n)),
    "Cnv": lambda spec, n: _twisted(spec, rotation.dihedral(n), _no_flip),
    "Cnh": _cnh,
    "S2n": _s2n,
    "Dn": lambda spec, n: _proper(spec, rotation.dihedral(n)),
    "Dnh": _dnh,
    "Dnd": _dnd,
}
POLYHEDRAL_FAMILY_TABLE: dict[str, Callable[[GroupSpec], PointGroup]] = {
    "T": lambda spec: _proper(spec, rotation.polyhedral("T")),
    "Td": lambda spec: _twisted(spec, rotation.polyhedral("O"), rotation.perm_even),
    "Th": lambda spec: _centro(spec, rotation.polyhedral("T")),
    "O": lambda spec: _proper(spec, rotation.polyhedral("O")),
    "Oh": lambda spec: _centro(spec, rotation.polyhedral("O")),
    "I": lambda spec: _proper(spec, rotation.polyhedral("I")),
    "Ih": lambda spec: _centro(spec, rotation.polyhedral("I")),
}


def point_group(spec: GroupSpec) -> PointGroup:
    """Build the abstract group of a parsed symbol by the family table."""
    if spec.n is None:
        return POLYHEDRAL_FAMILY_TABLE[spec.family](spec)
    return AXIAL_FAMILY_TABLE[spec.family](spec, spec.n)
