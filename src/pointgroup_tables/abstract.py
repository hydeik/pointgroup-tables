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


@dataclass(frozen=True)
class FamilyRule:
    """One row of the family table: for ``n`` of the given parity (``any``, ``even``
    or ``odd``), the group type and the rotation group ``R``, which is cyclic or
    dihedral of axis order ``n`` or, when ``doubled``, ``2n``. For ``S2n`` the
    ``n`` is that of the family name, half the subscript."""

    parity: str
    type: str
    kind: str  # "cyclic" or "dihedral"
    doubled: bool = False

    @property
    def rotation_group(self) -> str:
        """The name of ``R`` as the family table spells it: ``Cn``, ``C2n``, ``Dn``, ``D2n``."""
        return f"{'C' if self.kind == 'cyclic' else 'D'}{'2' if self.doubled else ''}n"

    def build(self, n: int) -> rotation.RotationGroup:
        order = 2 * n if self.doubled else n
        return (
            rotation.cyclic(order)
            if self.kind == "cyclic"
            else rotation.dihedral(order)
        )


@dataclass(frozen=True)
class PolyhedralRule:
    type: str
    rotation_group: str


# the family table of docs/spec/abstract-groups.md
AXIAL_FAMILY_TABLE: dict[str, tuple[FamilyRule, ...]] = {
    "Cn": (FamilyRule("any", PROPER, "cyclic"),),
    "Cnv": (FamilyRule("any", TWISTED, "dihedral"),),
    "Cnh": (
        FamilyRule("even", CENTRO, "cyclic"),
        FamilyRule("odd", TWISTED, "cyclic", doubled=True),
    ),
    "S2n": (
        FamilyRule("odd", CENTRO, "cyclic"),
        FamilyRule("even", TWISTED, "cyclic", doubled=True),
    ),
    "Dn": (FamilyRule("any", PROPER, "dihedral"),),
    "Dnh": (
        FamilyRule("even", CENTRO, "dihedral"),
        FamilyRule("odd", TWISTED, "dihedral", doubled=True),
    ),
    "Dnd": (
        FamilyRule("odd", CENTRO, "dihedral"),
        FamilyRule("even", TWISTED, "dihedral", doubled=True),
    ),
}
POLYHEDRAL_FAMILY_TABLE: dict[str, PolyhedralRule] = {
    "T": PolyhedralRule(PROPER, "T"),
    "Td": PolyhedralRule(TWISTED, "O"),
    "Th": PolyhedralRule(CENTRO, "T"),
    "O": PolyhedralRule(PROPER, "O"),
    "Oh": PolyhedralRule(CENTRO, "O"),
    "I": PolyhedralRule(PROPER, "I"),
    "Ih": PolyhedralRule(CENTRO, "I"),
}


def family_rule(family: str, n: int) -> FamilyRule:
    """The row of the family table that applies to ``n`` (the half subscript for ``S2n``)."""
    parity = "even" if n % 2 == 0 else "odd"
    for rule in AXIAL_FAMILY_TABLE[family]:
        if rule.parity in ("any", parity):
            return rule
    raise KeyError(family)  # pragma: no cover - every family covers both parities


def _build(
    spec: GroupSpec,
    group_type: str,
    R: rotation.RotationGroup,
    in_h: Callable[[Any], bool],
) -> PointGroup:
    """``G`` of the given type over ``R``; ``in_h`` is the index-2 subgroup of a twisted group."""
    if group_type == PROPER:
        return _proper(spec, R)
    if group_type == CENTRO:
        return _centro(spec, R)
    return _twisted(spec, R, in_h)


def _axial(spec: GroupSpec) -> PointGroup:
    """Build an axial group from its family rule. A twisted group's ``H`` is
    ``{k even}`` in a doubled ``R`` and ``{(k, 0)} = Cn`` in ``Dn`` (the ``Cnv`` case)."""
    n = spec.family_n
    rule = family_rule(spec.family, n)
    return _build(spec, rule.type, rule.build(n), _even_k if rule.doubled else _no_flip)


def _polyhedral(spec: GroupSpec) -> PointGroup:
    rule = POLYHEDRAL_FAMILY_TABLE[spec.family]
    # Td is the only twisted polyhedral group: H = the even permutations = T
    return _build(
        spec, rule.type, rotation.polyhedral(rule.rotation_group), rotation.perm_even
    )


def point_group(spec: GroupSpec) -> PointGroup:
    """Build the abstract group of a parsed symbol by the family table."""
    if spec.n is None:
        return _polyhedral(spec)
    return _axial(spec)
