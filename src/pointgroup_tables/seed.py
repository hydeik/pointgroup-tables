"""The seed table of known values (docs/spec/abstract-groups.md, docs/spec/double-groups.md).

One row per group of the ``seed`` set: the order, the class counts of ``G``
and ``G-dagger``, the number of projective irreps and the histogram of
element orders (group-theoretic order -> count). The rows were produced by
the prototypes under ``prototypes/`` and agree with the standard tables; the
``AG-2`` and ``DG-6`` checks compare a document with them.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class SeedRow:
    order: int
    class_count: int
    double_class_count: int
    projective_irreps: int
    element_orders: dict[int, int]


def _row(
    order: int,
    classes: int,
    double_classes: int,
    projective: int,
    orders: dict[str, int],
) -> SeedRow:
    return SeedRow(
        order,
        classes,
        double_classes,
        projective,
        {int(k): v for k, v in orders.items()},
    )


SEED: dict[str, SeedRow] = {
    "C1": _row(1, 1, 2, 1, {"1": 1}),
    "C1h": _row(2, 2, 4, 2, {"1": 1, "2": 1}),
    "S2": _row(2, 2, 4, 2, {"1": 1, "2": 1}),
    "C2": _row(2, 2, 4, 2, {"1": 1, "2": 1}),
    "C2v": _row(4, 4, 5, 1, {"1": 1, "2": 3}),
    "C2h": _row(4, 4, 8, 4, {"1": 1, "2": 3}),
    "D2": _row(4, 4, 5, 1, {"1": 1, "2": 3}),
    "D2h": _row(8, 8, 10, 2, {"1": 1, "2": 7}),
    "C3": _row(3, 3, 6, 3, {"1": 1, "3": 2}),
    "C3v": _row(6, 3, 6, 3, {"1": 1, "2": 3, "3": 2}),
    "C3h": _row(6, 6, 12, 6, {"1": 1, "2": 1, "3": 2, "6": 2}),
    "D3": _row(6, 3, 6, 3, {"1": 1, "2": 3, "3": 2}),
    "D3d": _row(12, 6, 12, 6, {"1": 1, "2": 7, "3": 2, "6": 2}),
    "D3h": _row(12, 6, 9, 3, {"1": 1, "2": 7, "3": 2, "6": 2}),
    "S6": _row(6, 6, 12, 6, {"1": 1, "2": 1, "3": 2, "6": 2}),
    "C4": _row(4, 4, 8, 4, {"1": 1, "2": 1, "4": 2}),
    "S4": _row(4, 4, 8, 4, {"1": 1, "2": 1, "4": 2}),
    "C4h": _row(8, 8, 16, 8, {"1": 1, "2": 3, "4": 4}),
    "C4v": _row(8, 5, 7, 2, {"1": 1, "2": 5, "4": 2}),
    "D4": _row(8, 5, 7, 2, {"1": 1, "2": 5, "4": 2}),
    "D2d": _row(8, 5, 7, 2, {"1": 1, "2": 5, "4": 2}),
    "D4h": _row(16, 10, 14, 4, {"1": 1, "2": 11, "4": 4}),
    "C6": _row(6, 6, 12, 6, {"1": 1, "2": 1, "3": 2, "6": 2}),
    "C6h": _row(12, 12, 24, 12, {"1": 1, "2": 3, "3": 2, "6": 6}),
    "C6v": _row(12, 6, 9, 3, {"1": 1, "2": 7, "3": 2, "6": 2}),
    "D6": _row(12, 6, 9, 3, {"1": 1, "2": 7, "3": 2, "6": 2}),
    "D6h": _row(24, 12, 18, 6, {"1": 1, "2": 15, "3": 2, "6": 6}),
    "T": _row(12, 4, 7, 3, {"1": 1, "2": 3, "3": 8}),
    "Td": _row(24, 5, 8, 3, {"1": 1, "2": 9, "3": 8, "4": 6}),
    "Th": _row(24, 8, 14, 6, {"1": 1, "2": 7, "3": 8, "6": 8}),
    "O": _row(24, 5, 8, 3, {"1": 1, "2": 9, "3": 8, "4": 6}),
    "Oh": _row(48, 10, 16, 6, {"1": 1, "2": 19, "3": 8, "4": 12, "6": 8}),
    "I": _row(60, 5, 9, 4, {"1": 1, "2": 15, "3": 20, "5": 24}),
    "Ih": _row(120, 10, 18, 8, {"1": 1, "2": 31, "3": 20, "5": 24, "6": 20, "10": 24}),
    "C5": _row(5, 5, 10, 5, {"1": 1, "5": 4}),
    "C5v": _row(10, 4, 8, 4, {"1": 1, "2": 5, "5": 4}),
    "D5h": _row(20, 8, 13, 5, {"1": 1, "2": 11, "5": 4, "10": 4}),
    "D5d": _row(20, 8, 16, 8, {"1": 1, "2": 11, "5": 4, "10": 4}),
    "C7v": _row(14, 5, 10, 5, {"1": 1, "2": 7, "7": 6}),
    "D7h": _row(28, 10, 17, 7, {"1": 1, "2": 15, "7": 6, "14": 6}),
    "D8d": _row(32, 11, 19, 8, {"1": 1, "2": 17, "4": 2, "8": 4, "16": 8}),
    "C8h": _row(16, 16, 32, 16, {"1": 1, "2": 3, "4": 4, "8": 8}),
    "S8": _row(8, 8, 16, 8, {"1": 1, "2": 1, "4": 2, "8": 4}),
    "D4d": _row(16, 7, 11, 4, {"1": 1, "2": 9, "4": 2, "8": 4}),
}
