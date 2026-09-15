"""Group sets and cached document models (docs/spec/verification.md, "Group sets")."""

from __future__ import annotations

import functools
from typing import Any

import pytest

from pointgroup_tables.document import model_document
from pointgroup_tables.seed import SEED
from pointgroup_tables.symbols import AXIAL_FAMILIES, axial_symbol, parse_group

LARGE_ORDERS = (13, 16, 17, 24, 30)


def _axial_symbols(orders: range | tuple[int, ...]) -> list[str]:
    symbols: list[str] = []
    for n in orders:
        for family in AXIAL_FAMILIES:
            subscript = 2 * n if family == "S2n" else n
            symbol = parse_group(axial_symbol(family, subscript)).symbol
            if symbol not in symbols:
                symbols.append(symbol)
    return symbols


SEED_GROUPS = list(SEED)
SWEEP_GROUPS = _axial_symbols(range(1, 13))
LARGE_GROUPS = _axial_symbols(LARGE_ORDERS)
DEFAULT_GROUPS = SEED_GROUPS + [symbol for symbol in SWEEP_GROUPS if symbol not in SEED]
POLYHEDRAL_GROUPS = ["T", "Td", "Th", "O", "Oh", "I", "Ih"]


@functools.cache
def model(symbol: str) -> dict[str, Any]:
    """The full document model of a group, built once per session."""
    return model_document(symbol)


def default_and_large() -> list[Any]:
    """``default`` plainly and ``large`` under the ``slow`` marker."""
    return [
        *DEFAULT_GROUPS,
        *(pytest.param(symbol, marks=pytest.mark.slow) for symbol in LARGE_GROUPS),
    ]
