"""The shipped Koster table: the Bethe symbols of the irreps of the 32
crystallographic point groups (docs/spec/irrep-labels.md, "Bethe aliases").

No rule generates Koster's ``Gamma_i`` numbering, so it is data: one entry
per irrep of each crystallographic group, carrying the group's identifying
elements (canonical labels), the irrep's characters on them in the value
grammar and its Bethe symbol. The file is generated from the transcriptions
under ``tests/reference/koster/bethe/`` by ``tools/koster_bethe.py
generate`` and never edited by hand. An entry is attached to the computed
irrep whose characters on the identifying elements it reproduces
(``REF-8``); a group outside the table receives no Bethe symbols.
"""

from __future__ import annotations

import tomllib
from functools import cache
from importlib import resources
from typing import Any

from .field import CyclotomicField, Value
from .irreps import Irrep
from .labels import Label
from .values import parse

KOSTER = "koster"


@cache
def koster_table() -> dict[str, Any]:
    """The parsed Koster table: ``name``, ``source`` and ``groups`` keyed by the
    normalised group symbol."""
    path = resources.files(__package__) / "data" / "conventions" / "bethe"
    return tomllib.loads((path / f"{KOSTER}.toml").read_text(encoding="utf-8"))


def bethe_groups() -> frozenset[str]:
    """The normalised symbols of the groups the Koster table covers."""
    return frozenset(koster_table()["groups"])


def attach_bethe(
    symbol: str, labels: list[Label], F: CyclotomicField, irreps: list[Irrep]
) -> None:
    """Set ``irrep.bethe`` on every irrep of a crystallographic group from the
    Koster table; leave the irreps of any other group untouched.

    Each entry names the irrep it belongs to by its characters on the group's
    identifying elements (and by ``projective``, which the characters on the
    unlifted elements of ``C1``, ``S2`` and ``C1h`` alone cannot separate).
    A table entry that matches no irrep, or two, or an irrep of another
    symbol, is corrupt shipped data and raises ``ValueError``.
    """
    entry = koster_table()["groups"].get(symbol)
    if entry is None:
        return
    if entry["conductor"] != F.conductor:
        raise ValueError(
            f"the Koster table of {symbol} lies in Q(zeta_{entry['conductor']}), "
            f"the document field is Q(zeta_{F.conductor})"
        )
    by_label = {label.text: k for k, label in enumerate(labels)}
    try:
        positions = [by_label[name] for name in entry["elements"]]
    except KeyError as exc:
        raise ValueError(
            f"the Koster table of {symbol} names an element {exc} the group has not"
        ) from None
    by_key: dict[tuple[bool, tuple[Value, ...]], Irrep] = {}
    for irrep in irreps:
        key = (irrep.projective, tuple(irrep.characters[k] for k in positions))
        if key in by_key:
            raise ValueError(
                f"the identifying elements of {symbol} do not separate "
                f"{by_key[key].symbol.text} from {irrep.symbol.text}"
            )
        by_key[key] = irrep
    for row in entry["irreps"]:
        key = (row["projective"], tuple(parse(F, text) for text in row["characters"]))
        irrep = by_key.get(key)
        if irrep is None:
            raise ValueError(
                f"the Koster table of {symbol} attaches {row['bethe']} to characters "
                f"{row['characters']} that no irrep has"
            )
        if irrep.symbol.text != row["symbol"]:
            raise ValueError(
                f"the Koster table of {symbol} names {row['symbol']} where the "
                f"characters give {irrep.symbol.text}"
            )
        if irrep.bethe is not None:
            raise ValueError(
                f"the Koster table of {symbol} attaches {irrep.bethe} and "
                f"{row['bethe']} to {irrep.symbol.text}"
            )
        irrep.bethe = row["bethe"]
    missing = [irrep.symbol.text for irrep in irreps if irrep.bethe is None]
    if missing:
        raise ValueError(
            f"the Koster table of {symbol} leaves {', '.join(missing)} without a Bethe symbol"
        )
