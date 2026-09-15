"""Conjugacy classes of ``G`` and their canonical labels (docs/spec/conjugacy-classes.md)."""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any

from .abstract import PointGroup
from .labels import Label, reduced_label


@dataclass(frozen=True)
class ConjugacyClass:
    index: int
    label: str
    members: tuple[int, ...]  # the representative first, the rest ascending
    inverse: int
    regular: bool

    @property
    def representative(self) -> int:
        return self.members[0]

    @property
    def size(self) -> int:
        return len(self.members)


def inverses(table: list[list[int]]) -> list[int]:
    """The inverse of every element, read off the Cayley table (identity is 0)."""
    return [row.index(0) for row in table]


def classes_from_table(table: list[list[int]]) -> list[list[int]]:
    """Classes of a group given by its Cayley table, ordered by first member."""
    n = len(table)
    inv = inverses(table)
    seen = [False] * n
    out = []
    for g in range(n):
        if seen[g]:
            continue
        members = sorted({table[table[h][g]][inv[h]] for h in range(n)})
        for x in members:
            seen[x] = True
        out.append(members)
    return out


def class_label(size: int, representative: Label, has_inverse: bool) -> str:
    return (str(size) if size > 1 else "") + reduced_label(representative, has_inverse)


def is_regular(G: PointGroup, element: tuple[Any, int]) -> bool:
    """Whether the class of ``element`` splits in the double group.

    The class of ``(r, p)`` splits exactly when the class of ``r`` splits in
    the binary cover of ``R``; a rotation's class merges with its barred
    element exactly when ``r`` is a binary rotation and ``R`` contains a binary
    rotation about a perpendicular axis, i.e. another element of order 2
    commuting with ``r``.
    """
    R = G.rotation_group
    r = element[0]
    if R.order_of(r) != 2:
        return True
    return not any(
        b != r and R.order_of(b) == 2 and R.mul(r, b) == R.mul(b, r) for b in R.elements
    )


def conjugacy_classes(
    G: PointGroup, table: list[list[int]], labels: list[Label]
) -> list[ConjugacyClass]:
    """The classes of ``G`` in canonical class order, labelled."""
    raw = classes_from_table(table)
    inv = inverses(table)
    class_of = {g: c for c, members in enumerate(raw) for g in members}
    out = []
    for c, members in enumerate(raw):
        representative = members[0]
        has_inverse = inv[representative] in members
        out.append(
            ConjugacyClass(
                index=c,
                label=class_label(len(members), labels[representative], has_inverse),
                members=tuple(members),
                inverse=class_of[inv[representative]],
                regular=is_regular(G, G.elements[representative]),
            )
        )
    return out


# ---- classes of the double group ------------------------------------------


@dataclass(frozen=True)
class DoubleClass:
    """A class of ``G-dagger``, ordered by base class (docs/spec/conjugacy-classes.md)."""

    index: int
    label: str
    members: tuple[int, ...]  # the representative first, the rest ascending
    inverse: int
    base_class: int
    split: bool
    barred: bool  # meaningful for split classes only

    @property
    def representative(self) -> int:
        return self.members[0]

    @property
    def size(self) -> int:
        return len(self.members)


def _with_representative_first(
    members: list[int], representative: int
) -> tuple[int, ...]:
    return (representative, *sorted(g for g in members if g != representative))


def double_class_label(
    size: int, representative: Label, has_inverse: bool, split: bool
) -> str:
    """A split class carries its full size; a merged one joins its two halves with ``+``."""
    if split:
        return class_label(size, representative, has_inverse)
    lifted = replace(representative, barred=False)
    barred = replace(representative, barred=True)
    half = size // 2
    return (
        class_label(half, lifted, has_inverse)
        + "+"
        + class_label(half, barred, has_inverse)
    )


def double_conjugacy_classes(
    table: list[list[int]],
    labels: list[Label],
    base_classes: list[ConjugacyClass],
) -> list[DoubleClass]:
    """The classes of ``G-dagger`` from its Cayley table, in the order by base class.

    ``labels`` are the labels of the ``2|G|`` elements of ``G-dagger`` (the barred
    block carries ``~``); the element ``i + |G|`` is the barred element of ``i``.
    """
    N = len(table) // 2
    raw = classes_from_table(table)
    inv = inverses(table)
    raw_of = {g: c for c, members in enumerate(raw) for g in members}
    ordered: list[
        tuple[list[int], int, int, bool, bool]
    ] = []  # members, rep, base, split, barred
    for base in base_classes:
        rep = base.representative
        lifted, barred = raw_of[rep], raw_of[rep + N]
        split = lifted != barred
        ordered.append((raw[lifted], rep, base.index, split, False))
        if split:
            ordered.append((raw[barred], rep + N, base.index, split, True))
    if len(ordered) != len(raw):
        raise ValueError(
            "the classes of the double group do not cover the base classes"
        )
    position = {raw_of[rep]: c for c, (_, rep, _, _, _) in enumerate(ordered)}
    out = []
    for c, (members, rep, base_index, split, barred) in enumerate(ordered):
        has_inverse = inv[rep] in members
        out.append(
            DoubleClass(
                index=c,
                label=double_class_label(len(members), labels[rep], has_inverse, split),
                members=_with_representative_first(members, rep),
                inverse=position[raw_of[inv[rep]]],
                base_class=base_index,
                split=split,
                barred=barred,
            )
        )
    return out
