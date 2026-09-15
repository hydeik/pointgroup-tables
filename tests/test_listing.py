"""The documents of the utility subcommands (docs/spec/cli.md, "Utility subcommands")."""

from __future__ import annotations

import json

import pytest
from test_elements import SCHEMA

from pointgroup_tables import abstract, conventions
from pointgroup_tables.document import Request, build_document
from pointgroup_tables.listing import (
    conventions_document,
    groups_document,
    shipped_schema,
)
from pointgroup_tables.symbols import ALIASES, AXIAL_KEYS, parse_group, subscript


def test_groups_document_layout() -> None:
    document = groups_document()
    assert list(document) == ["families", "polyhedral", "aliases"]
    assert [entry["key"] for entry in document["families"]] == list(AXIAL_KEYS)
    assert [entry["family"] for entry in document["families"]] == list(
        AXIAL_KEYS.values()
    )
    assert document["families"][1] == {
        "key": "Cv",
        "family": "Cnv",
        "constraint": "n >= 1",
        "rules": [{"n": "any", "type": "twisted", "rotation_group": "Dn"}],
        "example": "C3v",
    }
    assert document["families"][3] == {
        "key": "S",
        "family": "S2n",
        "constraint": "even subscript 2n >= 2",
        "rules": [
            {"n": "odd", "type": "centro", "rotation_group": "Cn"},
            {"n": "even", "type": "twisted", "rotation_group": "C2n"},
        ],
        "example": "S6",
    }
    assert document["polyhedral"] == [
        {"symbol": "T", "type": "proper", "rotation_group": "T"},
        {"symbol": "Td", "type": "twisted", "rotation_group": "O"},
        {"symbol": "Th", "type": "centro", "rotation_group": "T"},
        {"symbol": "O", "type": "proper", "rotation_group": "O"},
        {"symbol": "Oh", "type": "centro", "rotation_group": "O"},
        {"symbol": "I", "type": "proper", "rotation_group": "I"},
        {"symbol": "Ih", "type": "centro", "rotation_group": "I"},
    ]
    assert document["aliases"] == [
        {"input": alias, "symbol": symbol} for alias, symbol in ALIASES.items()
    ]
    assert [entry["input"] for entry in document["aliases"]] == [
        "Cs", "Ci", "C1v", "D1", "D1h", "D1d", "V", "Vh", "Vd",
    ]  # fmt: skip


@pytest.mark.parametrize("entry", groups_document()["families"], ids=lambda e: e["key"])
def test_family_rules_agree_with_the_constructed_groups(entry: dict) -> None:
    """Each rule of the listing is what building the group reports, for n = 1 .. 6."""
    for n in range(1, 7):
        symbol = f"{entry['key'][0]}{subscript(entry['family'], n)}{entry['key'][1:]}"
        spec = parse_group(symbol)
        if spec.family != entry["family"]:
            assert n == 1  # C1v, D1, D1h, D1d are aliases of other families
            continue
        rules = {rule["n"]: rule for rule in entry["rules"]}
        rule = rules.get("any") or rules["even" if n % 2 == 0 else "odd"]
        group = abstract.point_group(spec)
        assert group.type == rule["type"], (spec.symbol, rule)
        expected_order = {"Cn": n, "C2n": 2 * n, "Dn": 2 * n, "D2n": 4 * n}
        assert (
            len(group.rotation_group.elements) == expected_order[rule["rotation_group"]]
        )
        metadata = build_document(Request("elements", spec.symbol))["group"]
        assert metadata["family"] == entry["family"]
        assert metadata["rotation_group"] == group.rotation_group.name


def test_every_example_symbol_parses_to_its_family() -> None:
    for entry in groups_document()["families"]:
        assert parse_group(entry["example"]).family == entry["family"]


def test_conventions_document_lists_canonical_only_until_conventions_ship() -> None:
    document = conventions_document()
    assert list(document) == ["label_conventions", "basis_conventions"]
    assert document["label_conventions"] == [
        {"name": name, "source": block["source"]}
        for name, block in conventions.LABEL_CONVENTIONS.items()
    ]
    assert document["basis_conventions"][0] == {"name": "canonical"}
    assert all(
        list(entry) in (["name"], ["name", "source"])
        for entry in document["basis_conventions"]
    )


def test_shipped_schema_is_the_schema_file() -> None:
    data = shipped_schema()
    assert json.loads(data) == SCHEMA
    assert data.endswith(b"\n")
