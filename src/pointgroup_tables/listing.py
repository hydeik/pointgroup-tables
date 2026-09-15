"""The documents of the utility subcommands (docs/spec/cli.md, "Utility subcommands"):
``list groups``, ``list conventions`` and the shipped schema of ``schema``.

These documents are not table documents: ``schema.json`` does not cover them
and their layout is fixed by the command line chapter alone.
"""

from __future__ import annotations

from importlib import resources
from typing import Any

from . import conventions
from .abstract import AXIAL_FAMILY_TABLE, POLYHEDRAL_FAMILY_TABLE
from .symbols import ALIASES, AXIAL_KEYS, POLYHEDRAL_FAMILIES, axial_symbol, subscript

EXAMPLE_N = 3


def groups_document() -> dict[str, Any]:
    """What the group argument accepts: the axial families with their rules, the
    polyhedral groups and the aliases; no concrete axial group is enumerated."""
    families = []
    for key, family in AXIAL_KEYS.items():
        families.append(
            {
                "key": key,
                "family": family,
                "constraint": "even subscript 2n >= 2" if family == "S2n" else "n >= 1",
                "rules": [
                    {
                        "n": rule.parity,
                        "type": rule.type,
                        "rotation_group": rule.rotation_group,
                    }
                    for rule in AXIAL_FAMILY_TABLE[family]
                ],
                "example": axial_symbol(family, subscript(family, EXAMPLE_N)),
            }
        )
    polyhedral = [
        {
            "symbol": symbol,
            "type": POLYHEDRAL_FAMILY_TABLE[symbol].type,
            "rotation_group": POLYHEDRAL_FAMILY_TABLE[symbol].rotation_group,
        }
        for symbol in POLYHEDRAL_FAMILIES
    ]
    aliases = [{"input": alias, "symbol": symbol} for alias, symbol in ALIASES.items()]
    return {"families": families, "polyhedral": polyhedral, "aliases": aliases}


def _listed(table: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    """``{ name, source }`` per convention; the canonical basis has no source."""
    return [
        {key: block[key] for key in ("name", "source") if key in block}
        for block in table.values()
    ]


def conventions_document() -> dict[str, Any]:
    """What ``--label-convention`` and ``--basis-convention`` accept."""
    return {
        "label_conventions": _listed(conventions.LABEL_CONVENTIONS),
        "basis_conventions": _listed(conventions.BASIS_CONVENTIONS),
    }


def shipped_schema() -> bytes:
    """The bytes of the shipped ``schema.json``, as the ``schema`` subcommand writes them."""
    return (resources.files(__package__) / "data" / "schema.json").read_bytes()
