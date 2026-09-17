"""The ``altmann-herzig`` label convention applied to a document
(docs/spec/label-conventions.md, docs/spec/output-schema.md).

These checks use the shipped convention file only; the comparisons with the
transcriptions are the ``REF-*`` tests under the ``reference`` marker.
"""

from __future__ import annotations

import json

import jsonschema
import pytest
from test_elements import SCHEMA

from pointgroup_tables import conventions
from pointgroup_tables.cli import main
from pointgroup_tables.document import Request, build_document, model_document

AH = "altmann-herzig"


def test_the_shipped_label_convention_is_listed() -> None:
    assert list(conventions.LABEL_CONVENTIONS) == [AH]
    block = conventions.label_convention(AH)
    assert list(block) == ["name", "source"]
    assert block["name"] == AH
    assert "Altmann" in block["source"]


@pytest.mark.parametrize(
    ("symbol", "coverage"),
    [("C3v", "reference"), ("Oh", "reference"), ("D12h", "extrapolated")],
)
def test_coverage_follows_the_convention_range(symbol: str, coverage: str) -> None:
    document = build_document(Request("elements", symbol, label_convention=AH))
    assert list(document["label_convention"]) == ["name", "source", "coverage"]
    assert document["label_convention"]["coverage"] == coverage
    assert list(document)[:9] == [
        "schema_version", "generator", "group", "tables", "double", "values",
        "field", "label_convention", "basis_convention",
    ]  # fmt: skip


def test_the_key_is_absent_without_a_convention() -> None:
    assert "label_convention" not in model_document("C3v")
    assert "label_convention" in model_document("C3v", label_convention=AH)


def _labels(symbol: str, double: bool = False) -> dict[str, str]:
    """canonical label -> reference label of every element."""
    document = model_document(symbol, double=double, label_convention=AH)
    return {
        record["label"]: record["reference_label"]
        for record in document["elements"]["items"]
    }


def test_element_records_carry_the_reference_label_after_the_label() -> None:
    record = model_document("C3v", label_convention=AH)["elements"]["items"][1]
    assert list(record)[:4] == ["index", "label", "reference_label", "parity"]
    assert "reference_label" not in model_document("C3v")["elements"]["items"][1]


@pytest.mark.parametrize(
    ("symbol", "expected"),
    [
        # exponent_style = "plusminus" on the principal elements
        ("C8", {"C8^1": "C8+", "C8^-1": "C8-", "C8^3": "C8^3+", "C8^-3": "C8^3-"}),
        ("S10", {"S10^-3": "S10^3-", "S10^1": "S10+", "i": "i", "E": "E"}),
        # the master group's direction table: D4d inside D8h (spec example)
        ("D4d", {"C2'(1)": "C2'1", "C2'(2)": "C2'2", "C2''(1)": "C2'3", "C2''(2)": "C2'4"}),
        # the twist of C6v: sigma_v(1) keeps its letter (ADR-0007)
        ("C6v", {"sigma_v(1)": "sigma_v1", "sigma_d(1)": "sigma_d2"}),
        # per-group overrides win over every rule
        ("D2h", {"C2'(1)": "C2x", "sigma_v(1)": "sigma_y", "sigma_h": "sigma_z"}),
        ("C2v", {"sigma_v(1)": "sigma_y", "sigma_d(1)": "sigma_x", "C2": "C2"}),
        ("Th", {"S6^-1(1)": "S6,1-", "C3^1(3)": "C3,4+", "sigma_h(2)": "sigma_z"}),
        # a master beyond the table keeps the canonical indices (spelling rules only)
        ("D11", {"C2'(2)": "C2'(2)", "C11^-5": "C11^5-"}),
    ],
)  # fmt: skip
def test_reference_labels_of_elements(symbol: str, expected: dict[str, str]) -> None:
    labels = _labels(symbol)
    assert {label: labels[label] for label in expected} == expected


def swapped_elements(symbol: str) -> list[str]:
    """The canonical labels of ``G`` whose reference lift is our barred element:
    read off the double-group element list, where their unbarred lift carries
    the reference's tilde."""
    document = model_document(symbol, double=True, label_convention=AH)
    N = document["group"]["order"]
    return sorted(
        record["label"]
        for record in document["elements"]["items"][:N]
        if record["reference_label"].startswith("~")
    )


def test_double_group_labels_swap_the_tilde_of_the_listed_elements() -> None:
    labels = _labels("D3", double=True)
    assert labels["E"] == "E" and labels["~E"] == "~E"
    assert labels["C2'(1)"] == "C2'1" and labels["~C2'(1)"] == "~C2'1"
    # the book's C2'3 is the canonical C2'(2) (the [map] of the transcription)
    assert labels["C2'(2)"] == "~C2'3" and labels["~C2'(2)"] == "C2'3"
    assert labels["C3^-1"] == "C3-" and labels["~C3^-1"] == "~C3-"


@pytest.mark.parametrize(
    ("symbol", "expected"),
    [
        # docs/spec/double-groups.md, "Section overrides"
        ("D3", ["C2'(2)", "C2'(3)"]),
        ("O", ["C2'(1)", "C2'(3)", "C2'(5)"]),
        ("T", []),
        ("Th", []),
        # the lists tools/ah_transcription.py solves from the book's tilde marks
        # (tools/ah_convention.py check)
        ("C2v", ["sigma_d(1)"]),
        ("C3v", ["sigma_v(2)", "sigma_v(3)"]),
        ("D4", ["C2''(2)"]),
        ("D3h", ["C2'(2)", "C2'(3)", "sigma_v(3)"]),
        ("D4d", ["C2''(2)", "sigma_d(3)", "sigma_d(4)"]),
        ("C6v", ["sigma_d(2)", "sigma_d(3)", "sigma_v(2)", "sigma_v(3)"]),
        (
            "D6h",
            ["C2''(1)", "C2'(2)", "C2'(3)", "sigma_d(1)", "sigma_d(3)", "sigma_v(3)"],
        ),
        # no [secondary] entry for the master: nothing is swapped
        ("D11", []),
    ],
)
def test_swapped_lifts(symbol: str, expected: list[str]) -> None:
    assert swapped_elements(symbol) == expected


def test_canonical_fields_are_untouched_by_the_swaps() -> None:
    plain = model_document("D3", double=True)["elements"]
    labelled = model_document("D3", double=True, label_convention=AH)["elements"]

    def strip(record: dict) -> dict:
        return {k: v for k, v in record.items() if k != "reference_label"}

    assert [strip(r) for r in labelled["items"]] == plain["items"]
    assert [strip(r) for r in labelled["classes"]] == plain["classes"]


def _class_labels(symbol: str, double: bool = False) -> dict[str, str]:
    """canonical class label -> reference class label."""
    document = model_document(symbol, double=double, label_convention=AH)
    return {
        record["label"]: record["reference_label"]
        for record in document["elements"]["classes"]
    }


@pytest.mark.parametrize(
    ("symbol", "expected"),
    [
        # docs/spec/conjugacy-classes.md, "Reference class labels"
        ("Th", {"4C3": "4C3+", "4C3^-1": "4C3-", "4S6^-1": "4S6-", "4S6": "4S6+", "3sigma_h": "3sigma_h"}),
        ("D6h", {"2C3": "2C3", "2C6": "2C6", "3C2'": "3C2'", "3sigma_d": "3sigma_d"}),
        ("D2h", {"C2": "C2z", "C2'": "C2x", "C2''": "C2y", "sigma_v": "sigma_y", "sigma_d": "sigma_x", "sigma_h": "sigma_z"}),
        ("O", {"3C2": "3C2", "6C2'": "6C2'", "8C3": "8C3", "6C4": "6C4"}),
        ("Oh", {"3sigma_h": "3sigma_h", "6sigma_d": "6sigma_d", "8S6": "8S6", "6S4": "6S4"}),
        ("Ih", {"15sigma": "15sigma", "12S10^3": "12S10^3", "12C5^2": "12C5^2"}),
        ("S8", {"S8": "S8+", "S8^-1": "S8-", "S8^3": "S8^3+", "S8^-3": "S8^3-"}),
    ],
)  # fmt: skip
def test_reference_class_labels(symbol: str, expected: dict[str, str]) -> None:
    labels = _class_labels(symbol)
    assert {label: labels[label] for label in expected} == expected


def test_reference_class_labels_of_the_double_group() -> None:
    labels = _class_labels("D3", double=True)
    assert labels == {
        "E": "E", "~E": "~E", "2C3": "2C3", "2~C3": "2~C3",
        "3C2'": "3C2'", "3~C2'": "3~C2'",
    }  # fmt: skip
    assert _class_labels("T", double=True)["4~C3^-1"] == "4~C3-"
    # merged classes: the reference's unbarred half first
    assert _class_labels("C2v", double=True)["sigma_d+~sigma_d"] == "sigma_x+~sigma_x"
    assert _class_labels("D4", double=True)["2C2''+2~C2''"] == "2C2''+2~C2''"


def test_class_records_carry_the_reference_label_after_the_label() -> None:
    document = model_document("C3v", label_convention=AH)
    record = document["elements"]["classes"][1]
    assert list(record)[:4] == ["index", "label", "reference_label", "size"]
    # the character table's columns are the classes of G-dagger, labelled too
    columns = document["characters"]["classes"]
    assert [column["reference_label"] for column in columns] == [
        "E", "~E", "2C3", "2~C3", "3sigma_v", "3~sigma_v"
    ]  # fmt: skip


def _irrep_symbols(symbol: str) -> dict[str, str]:
    document = model_document(symbol, label_convention=AH)
    return {
        record["symbol"]: record["reference_symbol"]
        for record in document["characters"]["irreps"]
    }


@pytest.mark.parametrize(
    ("symbol", "expected"),
    [
        # docs/spec/irrep-labels.md, "Reference irrep labels in a label convention"
        ("C2h", {"Ag": "Ag", "1E1/2g": "1E1/2,g", "2E1/2u": "2E1/2,u"}),
        ("S6", {"B3/2g": "A3/2,g", "B3/2u": "A3/2,u"}),
        ("C3", {"B3/2": "A3/2", "1E1/2": "1E1/2"}),
        # the book's member 1 of a complex pair is the canonical member 2 in the
        # cyclic-type groups (decided by the characters of T 3.4, T 13.4, ...)
        ("C3", {"1E": "2E", "2E": "1E"}),
        ("S6", {"1Eg": "2Eg", "2Eu": "1Eu"}),
        ("C4", {"1E1/2": "2E1/2", "2E3/2": "2E3/2"}),
        # the book indexes the E_k of D4d differently; the E_k/2 of D5h agree once the
        # interchanged symbols of T 34.4 are corrected (tests/reference/altmann-herzig/
        # characters/D5h.toml, decided by T 34.6 and T 34.7)
        ("D4d", {"E1": "E3", "E3": "E1", "E2": "E2"}),
        ("D5h", {"E3/2": "E3/2", "E7/2": "E7/2", "E1/2": "E1/2"}),
        ("T", {"1E3/2": "1F3/2", "2E3/2": "2F3/2", "E1/2": "E1/2", "T": "T"}),
        ("Th", {"1E3/2g": "1F3/2,g", "2E3/2u": "2F3/2,u", "Tg": "Tg"}),
        ("O", {"G3/2": "F3/2", "T2": "T2"}),
        ("I", {"G": "F", "G3/2": "F3/2", "I5/2": "I5/2", "H": "H"}),
        ("Ih", {"Gu": "Fu", "G3/2g": "F3/2,g", "I5/2u": "I5/2,u"}),
        ("C1h", {"A'": "A'", "A''": "A''", "1E1/2": "1E1/2"}),
        ("D3h", {"E''": "E''", "E5/2": "E5/2"}),
        ("D2d", {"B1": "B1", "E3/2": "E3/2"}),
    ],
)  # fmt: skip
def test_reference_irrep_symbols(symbol: str, expected: dict[str, str]) -> None:
    symbols = _irrep_symbols(symbol)
    assert {symbol: symbols[symbol] for symbol in expected} == expected


def test_irrep_records_carry_the_reference_symbol_after_the_symbol() -> None:
    document = model_document("O", label_convention=AH)
    for table in ("characters", "matrices"):
        record = document[table]["irreps"][7]
        assert list(record)[:4] == ["index", "symbol", "reference_symbol", "bethe"]
        assert record["reference_symbol"] == "F3/2"


def test_documents_with_a_label_convention_validate_against_the_schema() -> None:
    for symbol in ("C3v", "D2h", "Oh"):
        for subcommand, double in (("elements", True), ("all", False)):
            request = Request(subcommand, symbol, label_convention=AH, double=double)
            jsonschema.Draft202012Validator(SCHEMA).validate(build_document(request))


def _run(
    capsysbinary: pytest.CaptureFixture[bytes], *argv: str
) -> tuple[int, bytes, str]:
    code = main(list(argv))
    captured = capsysbinary.readouterr()
    return code, captured.out, captured.err.decode()


def test_cli_rejects_an_unknown_label_convention(
    capsysbinary: pytest.CaptureFixture[bytes],
) -> None:
    code, out, err = _run(capsysbinary, "elements", "C3v", "--label-convention", "nope")
    assert code == 4
    assert out == b""
    assert err == (
        "pointgroup-tables: error: unknown label convention 'nope'\n"
        "hint: run 'pointgroup-tables list conventions' for the shipped conventions\n"
    )


def test_cli_lists_and_applies_the_shipped_label_convention(
    capsysbinary: pytest.CaptureFixture[bytes],
) -> None:
    code, out, err = _run(capsysbinary, "list", "conventions")
    assert (code, err) == (0, "")
    listed = json.loads(out)["label_conventions"]
    assert [entry["name"] for entry in listed] == [AH]
    code, out, err = _run(
        capsysbinary, "elements", "Oh", "--double-group", "--label-convention", AH
    )
    assert (code, err) == (0, "")
    document = json.loads(out)
    assert document["label_convention"] == {
        "name": AH,
        "source": listed[0]["source"],
        "coverage": "reference",
    }
    labels = {r["label"]: r["reference_label"] for r in document["elements"]["items"]}
    assert labels["C2'(1)"] == "~C2'd" and labels["~C2'(1)"] == "C2'd"
    assert labels["C2'(2)"] == "C2'c" and labels["~C4^1(1)"] == "~C4z+"
