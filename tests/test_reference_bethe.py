"""``REF-8`` and the Bethe aliases of the irrep record (docs/spec/irrep-labels.md,
"Bethe aliases"; docs/spec/verification.md, "Reference comparison").

The shipped Koster table attaches a Bethe symbol to every irrep of the 32
crystallographic groups by the irrep's characters on the group's identifying
elements. ``REF-8`` reads the transcriptions under
``tests/reference/koster/bethe/`` and checks, through the document, that the
irrep carrying each ``Gamma_i`` has the transcribed characters. The other
tests cover the shipped table itself: ``bethe`` is present for exactly the 32
groups, distinct within a group, grammar-conformant, and carries the parity
sign exactly in the centrosymmetric groups, agreeing with the ``g``/``u``
suffix.
"""

from __future__ import annotations

import re
import tomllib
from fractions import Fraction
from pathlib import Path
from typing import Any

import pytest
from conftest import DEFAULT_GROUPS, model

from pointgroup_tables.document import Request, build_document
from pointgroup_tables.field import CyclotomicField, Value, cyclotomic_field
from pointgroup_tables.values import parse

KOSTER = Path(__file__).resolve().parent / "reference" / "koster" / "bethe"
FILES = sorted(KOSTER.glob("*.toml")) if KOSTER.is_dir() else []

# the 32 crystallographic groups by normalised symbol (Cs = C1h, Ci = S2)
CRYSTALLOGRAPHIC = frozenset(
    {
        "C1", "S2", "C2", "C1h", "C2h", "D2", "C2v", "D2h",
        "C4", "S4", "C4h", "D4", "C4v", "D2d", "D4h",
        "C3", "S6", "D3", "C3v", "D3d",
        "C6", "C3h", "C6h", "D6", "C6v", "D3h", "D6h",
        "T", "Th", "O", "Td", "Oh",
    }
)  # fmt: skip
BETHE = re.compile(r"^Gamma([1-9][0-9]*)([+-])?$")

# the transcription templates name a root of unity by its own order (zeta3^2);
# the value grammar names it by the document conductor (zeta^4 in Q(zeta_6))
_TRANSCRIBED_ROOT = re.compile(r"^(-)?(?:(\d+)\*)?zeta(\d+)(?:\^(\d+))?$")


def _transcribed(F: CyclotomicField, text: str) -> Value:
    match = _TRANSCRIBED_ROOT.match(text)
    if match is None:
        return parse(F, text)
    sign, factor, order, power = match.groups()
    value = F.root_of_unity(int(order), int(power or 1))
    if factor:
        value = F.scale(value, Fraction(int(factor)))
    return F.neg(value) if sign else value


def _column_of_element(document: dict[str, Any], label: str) -> int:
    """The column of the character table (classes of ``G-dagger``) holding the
    unbarred lift of the element with this canonical label."""
    items = document["elements"]["items"]
    [item] = [item for item in items if item["label"] == label]
    columns = document["characters"]["classes"]
    [column] = [
        c
        for c, cls in enumerate(columns)
        if cls["base_class"] == item["class"] and not cls.get("barred", False)
    ]
    return column


# ---- REF-8 -------------------------------------------------------------------------


@pytest.mark.reference
@pytest.mark.parametrize("path", FILES, ids=lambda path: path.stem)
def test_ref8_each_bethe_symbol_sits_on_the_irrep_with_its_characters(
    path: Path,
) -> None:
    transcription = tomllib.loads(path.read_text())
    document = model(transcription["group"])
    assert document["group"]["symbol"] == path.stem
    F = cyclotomic_field(document["field"]["conductor"])
    columns = [_column_of_element(document, name) for name in transcription["elements"]]
    table = document["characters"]
    by_bethe: dict[str, dict[str, Any]] = {}
    for record in table["irreps"]:
        assert record["bethe"] not in by_bethe
        by_bethe[record["bethe"]] = record
    assert len(transcription["irreps"]) == len(table["irreps"])
    for entry in transcription["irreps"]:
        record = by_bethe[entry["bethe"]]
        row = table["exact"][record["index"]]
        expected = [F.exact(_transcribed(F, text)) for text in entry["characters"]]
        assert [row[c] for c in columns] == expected, entry["symbol"]
        assert record["projective"] == entry["projective"]
        assert record["dimension"] == entry["dimension"]
        assert record["symbol"] == entry["symbol"]


@pytest.mark.reference
def test_every_crystallographic_group_is_transcribed() -> None:
    if not FILES:
        pytest.skip(f"no transcription under {KOSTER}")
    assert {path.stem for path in FILES} == CRYSTALLOGRAPHIC


# ---- the shipped table ---------------------------------------------------------------


@pytest.mark.parametrize("symbol", DEFAULT_GROUPS)
def test_bethe_is_present_for_exactly_the_crystallographic_groups(symbol: str) -> None:
    records = model(symbol)["characters"]["irreps"]
    present = {"bethe" in record for record in records}
    assert present == {symbol in CRYSTALLOGRAPHIC}


@pytest.mark.parametrize("symbol", sorted(CRYSTALLOGRAPHIC))
def test_bethe_symbols_are_distinct_and_signed_by_parity(symbol: str) -> None:
    document = model(symbol)
    centrosymmetric = document["group"]["type"] == "centro"
    records = document["characters"]["irreps"]
    symbols = [record["bethe"] for record in records]
    assert len(set(symbols)) == len(symbols)
    for record in records:
        match = BETHE.match(record["bethe"])
        assert match is not None, record["bethe"]
        sign = match.group(2)
        assert (sign is not None) == centrosymmetric
        if sign is not None:
            assert (sign == "+") == (record["suffix"] == "g")


def test_bethe_is_emitted_by_the_characters_subcommand_and_survives_irreps() -> None:
    # the examples of the irrep-labels chapter: Koster's numbers do not follow
    # the canonical order (D3h: E' is Gamma6, E'' is Gamma5), and the parity
    # sign of a centrosymmetric group rides on the number (Oh, D3d)
    examples = {
        ("D3h", "all"): {"E'": "Gamma6", "E''": "Gamma5", "E5/2": "Gamma8"},
        ("Oh", "ordinary"): {"T1g": "Gamma4+", "T2u": "Gamma5-"},
        ("D3d", "projective"): {"1E3/2g": "Gamma5+", "2E3/2u": "Gamma6-"},
        ("C6v", "all"): {"B1": "Gamma4", "B2": "Gamma3", "E1": "Gamma5"},
    }
    for (group, which), expected in examples.items():
        document = build_document(Request("characters", group, irreps=which))
        bethe = {r["symbol"]: r["bethe"] for r in document["characters"]["irreps"]}
        for symbol, gamma in expected.items():
            assert bethe[symbol] == gamma, (group, symbol)


def test_bethe_sits_between_symbol_and_dimension_in_the_record() -> None:
    document = build_document(Request("characters", "Ci"))
    for record in document["characters"]["irreps"]:
        keys = list(record)
        assert keys[:4] == ["index", "symbol", "bethe", "dimension"]
