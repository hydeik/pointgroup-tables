"""The ``elements`` document: group and class checks, labels, schema (``AG-*``, ``CL-*``, ``OS-1``)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import jsonschema
import pytest
from conftest import (
    DEFAULT_GROUPS,
    LARGE_GROUPS,
    POLYHEDRAL_GROUPS,
    SEED_GROUPS,
    default_and_large,
    model,
)

from pointgroup_tables import checks
from pointgroup_tables.document import Request, build_document, parse_values
from pointgroup_tables.seed import SEED

REPO = Path(__file__).resolve().parent.parent
SCHEMA = json.loads(
    (REPO / "src" / "pointgroup_tables" / "data" / "schema.json").read_text()
)
SAMPLES = REPO / "prototypes" / "samples"


def labels_of(symbol: str) -> list[str]:
    return [record["label"] for record in model(symbol)["elements"]["items"]]


def class_labels_of(symbol: str) -> list[str]:
    return [record["label"] for record in model(symbol)["elements"]["classes"]]


# ---- AG-* and CL-* as pure functions of the model --------------------------


@pytest.mark.parametrize(
    "check_id", ["AG-1", "AG-3", "CL-1", "CL-2", "CL-3", "CL-4", "CL-5"]
)
@pytest.mark.parametrize("symbol", default_and_large())
def test_group_and_class_checks(check_id: str, symbol: str) -> None:
    assert checks.check(model(symbol), [check_id]) == []


@pytest.mark.parametrize("symbol", SEED_GROUPS)
def test_ag2_seed_table(symbol: str) -> None:
    assert checks.check(model(symbol), ["AG-2"]) == []


@pytest.mark.parametrize("symbol", SEED_GROUPS)
def test_regular_classes_count_projective_irreps(symbol: str) -> None:
    regular = sum(1 for cls in model(symbol)["elements"]["classes"] if cls["regular"])
    assert regular == SEED[symbol].projective_irreps


def test_seed_table_covers_the_seed_set() -> None:
    crystallographic = {
        "C1", "S2", "C2", "C1h", "C2h", "D2", "C2v", "D2h", "C4", "S4", "C4h", "D4", "C4v", "D2d",
        "D4h", "C3", "S6", "D3", "C3v", "D3d", "C6", "C3h", "C6h", "D6", "C6v", "D3h", "D6h",
        "T", "Th", "O", "Td", "Oh",
    }  # fmt: skip
    extra = {
        "I",
        "Ih",
        "C5",
        "C5v",
        "D5h",
        "D5d",
        "C7v",
        "D7h",
        "D8d",
        "C8h",
        "S8",
        "D4d",
    }
    assert set(SEED) == crystallographic | extra


def test_checks_refuse_a_document_without_the_needed_table() -> None:
    document = build_document(Request("elements", "C3v"))
    assert "cayley" not in document
    assert checks.applicable_checks(document) == [
        "AG-3", "CL-2", "CL-4", "CL-5", "DG-6",
    ]  # fmt: skip
    with pytest.raises(ValueError, match="AG-1 needs the cayley table"):
        checks.check(document, ["AG-1"])


# ---- canonical labels and class labels -------------------------------------


@pytest.mark.parametrize(
    ("symbol", "expected"),
    [
        ("C3v", "E C3^1 C3^-1 sigma_v(1) sigma_v(2) sigma_v(3)"),
        (
            "D6h",
            (
                "E C6^1 C3^1 C2 C3^-1 C6^-1 C2'(1) C2''(1) C2'(2) C2''(2) C2'(3) C2''(3) "
                "i S3^-1 S6^-1 sigma_h S6^1 S3^1 sigma_d(2) sigma_v(3) sigma_d(3) sigma_v(1) sigma_d(1) sigma_v(2)"
            ),
        ),
        ("S8", "E C4^1 C2 C4^-1 S8^-3 S8^-1 S8^1 S8^3"),
        (
            "D5d",
            (
                "E C5^1 C5^2 C5^-2 C5^-1 C2'(1) C2'(2) C2'(3) C2'(4) C2'(5) "
                "i S10^-3 S10^-1 S10^1 S10^3 sigma_d(1) sigma_d(2) sigma_d(3) sigma_d(4) sigma_d(5)"
            ),
        ),
        ("C6", "E C6^1 C3^1 C2 C3^-1 C6^-1"),
        ("D2d", "E C2 C2'(1) C2''(1) S4^-1 S4^1 sigma_d(1) sigma_d(2)"),
        ("C3h", "E C3^1 C3^-1 S3^-1 sigma_h S3^1"),
        ("C1h", "E sigma_h"),
        ("S2", "E i"),
    ],
)
def test_canonical_element_labels(symbol: str, expected: str) -> None:
    assert labels_of(symbol) == expected.split()


@pytest.mark.parametrize(
    ("symbol", "expected"),
    [
        ("C3v", "E 2C3 3sigma_v"),
        ("D6h", "E 2C6 2C3 C2 3C2' 3C2'' i 2S3 2S6 sigma_h 3sigma_d 3sigma_v"),
        ("D3h", "E 2C3 3C2' 2S3 sigma_h 3sigma_v"),
        ("Oh", "E 3C2 6C2' 8C3 6C4 i 3sigma_h 6sigma_d 8S6 6S4"),
        ("Ih", "E 15C2 20C3 12C5 12C5^2 i 15sigma 20S6 12S10^3 12S10"),
        ("Td", "E 3C2 8C3 6sigma_d 6S4"),
        ("Th", "E 3C2 4C3 4C3^-1 i 3sigma_h 4S6^-1 4S6"),
        ("D2d", "E C2 2C2' 2S4 2sigma_d"),
        ("S8", "E C4 C2 C4^-1 S8^-3 S8^-1 S8 S8^3"),
        ("D5d", "E 2C5 2C5^2 5C2' i 2S10^3 2S10 5sigma_d"),
        ("C6", "E C6 C3 C2 C3^-1 C6^-1"),
    ],
)
def test_canonical_class_labels(symbol: str, expected: str) -> None:
    assert class_labels_of(symbol) == expected.split()


@pytest.mark.parametrize("symbol", default_and_large())
def test_element_labels_match_the_grammar_and_are_distinct(symbol: str) -> None:
    labels = labels_of(symbol)
    assert all(checks.ELEMENT_LABEL.match(label) for label in labels)
    assert len(set(labels)) == len(labels)


def test_ih_polyhedral_records() -> None:
    items = model("Ih")["elements"]["items"]
    c5 = next(record for record in items if record["label"] == "C5^1(1)")
    assert c5["kind"] == "rotation" and c5["order"] == 5 and c5["power"] == 1
    assert c5["axis_class"] == "polyhedral" and c5["axis_index"] == 1
    assert items[c5["inverse"]]["label"] == "C5^-1(1)"
    s10 = next(record for record in items if record["label"] == "S10^-3(1)")
    assert s10["kind"] == "rotoreflection" and s10["order"] == 10 and s10["power"] == -3
    sigma = next(record for record in items if record["label"] == "sigma(1)")
    assert (
        sigma["kind"] == "reflection" and sigma["order"] == 2 and "power" not in sigma
    )
    inversion = next(record for record in items if record["label"] == "i")
    assert inversion["kind"] == "inversion" and "axis_class" not in inversion


# ---- the prototype samples as a regression oracle --------------------------


@pytest.mark.parametrize(
    ("symbol", "sample"),
    [
        ("C3v", "C3v_elements.json"),
        ("D6h", "D6h_all.json"),
        ("Oh", "Oh_all.json"),
        ("Ih", "Ih_elements.json"),
    ],
)
def test_elements_table_equals_the_prototype_sample(symbol: str, sample: str) -> None:
    expected = json.loads((SAMPLES / sample).read_text())
    document = build_document(Request("elements", symbol))
    if symbol in POLYHEDRAL_GROUPS:
        # the prototype read the axis class off the label's spelling; every axis of
        # a polyhedral group is a polyhedral axis (docs/spec/label-conventions.md, rule 5)
        for record in expected["elements"]["items"]:
            if "axis_class" in record:
                record["axis_class"] = "polyhedral"
    assert document["elements"] == expected["elements"]
    assert document["group"] == expected["group"]
    assert document["field"] == expected["field"]


# ---- the document: metadata and schema (OS-1) ------------------------------


def test_metadata_keys_in_fixed_order() -> None:
    document = build_document(Request("elements", "c3v"))
    assert list(document) == [
        "schema_version", "generator", "group", "tables", "double", "values", "field", "basis_convention", "elements",
    ]  # fmt: skip
    assert document["group"] == {
        "symbol": "C3v", "input_symbol": "c3v", "family": "Cnv", "n": 3,
        "order": 6, "type": "twisted", "rotation_group": "D3",
    }  # fmt: skip
    assert document["tables"] == ["elements"]
    assert document["double"] is False
    assert document["values"] == ["exact", "numeric"]
    assert document["basis_convention"] == {"name": "canonical"}
    assert "label_convention" not in document
    assert "irreps" not in document


def test_polyhedral_metadata_has_no_n() -> None:
    document = build_document(Request("elements", "Oh"))
    assert document["group"] == {
        "symbol": "Oh", "input_symbol": "Oh", "family": "Oh", "order": 48, "type": "centro", "rotation_group": "O",
    }  # fmt: skip
    assert document["field"]["conductor"] == 8


@pytest.mark.parametrize(
    ("text", "forms"),
    [
        ("both", ("exact", "numeric")),
        ("exact", ("exact",)),
        ("pretty,numeric", ("numeric", "pretty")),
        ("both,pretty", ("exact", "numeric", "pretty")),
    ],
)
def test_parse_values(text: str, forms: tuple[str, ...]) -> None:
    assert parse_values(text) == forms
    assert build_document(Request("elements", "C1", values=forms))["values"] == list(
        forms
    )


VALUE_SETS = [
    "exact",
    "numeric",
    "both",
    "pretty",
    "exact,pretty",
    "numeric,pretty",
    "both,pretty",
]


@pytest.mark.parametrize("values", VALUE_SETS)
@pytest.mark.parametrize("symbol", SEED_GROUPS)
def test_elements_document_validates_against_the_schema(
    symbol: str, values: str
) -> None:
    document = build_document(Request("elements", symbol, values=parse_values(values)))
    jsonschema.Draft202012Validator(SCHEMA).validate(document)


@pytest.mark.parametrize(
    "symbol",
    [
        *(symbol for symbol in DEFAULT_GROUPS if symbol not in SEED),
        *(pytest.param(symbol, marks=pytest.mark.slow) for symbol in LARGE_GROUPS),
    ],
)
def test_sweep_and_large_documents_validate_against_the_schema(symbol: str) -> None:
    jsonschema.Draft202012Validator(SCHEMA).validate(
        build_document(Request("elements", symbol))
    )


def test_shipped_schema_equals_the_specification_schema() -> None:
    spec = json.loads((REPO / "docs" / "spec" / "schema.json").read_text())
    assert SCHEMA == spec


# ---- class records ------------------------------------------------------------


@pytest.mark.parametrize("symbol", POLYHEDRAL_GROUPS + ["C3v", "D6h", "C4h"])
def test_class_records_are_self_consistent(symbol: str) -> None:
    document: dict[str, Any] = model(symbol)
    items = document["elements"]["items"]
    for cls in document["elements"]["classes"]:
        assert cls["members"][0] == cls["representative"]
        assert cls["members"][1:] == sorted(cls["members"][1:])
        assert all(items[g]["class"] == cls["index"] for g in cls["members"])
        assert list(cls) == [
            "index",
            "label",
            "size",
            "representative",
            "members",
            "inverse",
            "regular",
        ]
