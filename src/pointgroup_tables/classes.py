"""Conjugacy classes of ``G`` and their canonical labels (docs/spec/conjugacy-classes.md)."""

from __future__ import annotations

from dataclasses import dataclass
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
