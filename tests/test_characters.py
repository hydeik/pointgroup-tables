"""The ``characters`` document of the axial families (docs/spec/irrep-labels.md,
docs/spec/conjugacy-classes.md, docs/spec/output-schema.md): ``IL-*``, ``CH-*``,
``VL-*``, the symbol rules, ``--irreps``, ``--values`` and the schema."""

from __future__ import annotations

import itertools
import json
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
from test_elements import SAMPLES, SCHEMA

from pointgroup_tables import checks
from pointgroup_tables.document import (
    VALUE_FORMS,
    Request,
    build_document,
    build_model,
    model_document,
    parse_values,
)
from pointgroup_tables.errors import PointGroupTablesError, UsageError
from pointgroup_tables.field import cyclotomic_field
from pointgroup_tables.irreps import Symbol
from pointgroup_tables.symbols import parse_group

AXIAL_SEED = [symbol for symbol in SEED_GROUPS if symbol not in POLYHEDRAL_GROUPS]
CHECKS = [
    "IL-1",
    "IL-2",
    "IL-3",
    "IL-4",
    "CH-1",
    "CH-2",
    "CH-3",
    "CH-4",
    "CH-5",
    "VL-1",
    "VL-2",
]


def characters_of(symbol: str) -> dict[str, Any]:
    document = model(symbol)
    if "characters" not in document:
        pytest.skip(
            f"the character table of {symbol} is not built yet (polyhedral families)"
        )
    return document["characters"]


def symbols_of(symbol: str) -> tuple[str, str]:
    table = characters_of(symbol)
    ordinary = " ".join(r["symbol"] for r in table["irreps"] if not r["projective"])
    projective = " ".join(r["symbol"] for r in table["irreps"] if r["projective"])
    return ordinary, projective


# ---- the ledger checks as pure functions of the model ----------------------------------


@pytest.mark.parametrize("check_id", CHECKS)
@pytest.mark.parametrize("symbol", default_and_large())
def test_irrep_and_character_checks(check_id: str, symbol: str) -> None:
    characters_of(symbol)
    assert checks.check(model(symbol), [check_id]) == []


def _ch6_groups() -> list[Any]:
    """``CH-6`` runs plainly on ``|G-dagger| <= 96`` and under ``slow`` beyond."""
    out: list[Any] = []
    for symbol in DEFAULT_GROUPS + LARGE_GROUPS:
        spec = parse_group(symbol)
        if spec.polyhedral:
            continue
        order = build_model(spec).order
        if 2 * order > 96 or symbol in LARGE_GROUPS:
            out.append(pytest.param(symbol, marks=pytest.mark.slow))
        else:
            out.append(symbol)
    return out


@pytest.mark.parametrize("symbol", _ch6_groups())
def test_ch6_class_algebra_homomorphisms(symbol: str) -> None:
    assert checks.check(model(symbol), ["CH-6"]) == []


@pytest.mark.parametrize("which", ["ordinary", "projective"])
@pytest.mark.parametrize("symbol", AXIAL_SEED)
def test_checks_hold_on_filtered_documents(symbol: str, which: str) -> None:
    """Row and column orthogonality, the block sums and the class-function property
    hold on the filtered rows with their partial expectations."""
    document = build_document(Request("characters", symbol, irreps=which))
    applicable = [
        name
        for name in checks.applicable_checks(document)
        if name.startswith(("IL", "CH", "VL"))
    ]
    assert {
        "IL-1",
        "IL-2",
        "IL-4",
        "CH-1",
        "CH-2",
        "CH-3",
        "CH-4",
        "CH-5",
        "VL-1",
        "VL-2",
    } <= set(applicable)
    assert checks.check(document, applicable) == []


# ---- symbols and order (docs/spec/irrep-labels.md, "Examples") ---------------------------

EXAMPLES = [
    ("C1", "A", "A1/2"),
    ("Cs", "A' A''", "1E1/2 2E1/2"),
    ("Ci", "Ag Au", "A1/2g A1/2u"),
    ("C2", "A B", "1E1/2 2E1/2"),
    ("C3", "A 1E 2E", "B3/2 1E1/2 2E1/2"),
    ("C4", "A B 1E 2E", "1E1/2 2E1/2 1E3/2 2E3/2"),
    ("S4", "A B 1E 2E", "1E1/2 2E1/2 1E3/2 2E3/2"),
    ("S6", "Ag 1Eg 2Eg Au 1Eu 2Eu", "B3/2g 1E1/2g 2E1/2g B3/2u 1E1/2u 2E1/2u"),
    ("C3h", "A' 1E' 2E' A'' 1E'' 2E''", "1E1/2 2E1/2 1E3/2 2E3/2 1E5/2 2E5/2"),
    ("C5", "A 1E1 2E1 1E2 2E2", "B5/2 1E1/2 2E1/2 1E3/2 2E3/2"),
    ("C6", "A B 1E1 2E1 1E2 2E2", "1E1/2 2E1/2 1E3/2 2E3/2 1E5/2 2E5/2"),
    ("D2", "A B1 B2 B3", "E1/2"),
    ("D2h", "Ag B1g B2g B3g Au B1u B2u B3u", "E1/2g E1/2u"),
    ("C2v", "A1 A2 B1 B2", "E1/2"),
    ("C3v", "A1 A2 E", "E1/2 1E3/2 2E3/2"),
    ("D3", "A1 A2 E", "E1/2 1E3/2 2E3/2"),
    ("C4v", "A1 A2 B1 B2 E", "E1/2 E3/2"),
    ("D4", "A1 A2 B1 B2 E", "E1/2 E3/2"),
    ("D2d", "A1 A2 B1 B2 E", "E1/2 E3/2"),
    ("D4h", "A1g A2g B1g B2g Eg A1u A2u B1u B2u Eu", "E1/2g E3/2g E1/2u E3/2u"),
    ("D3h", "A1' A2' E' A1'' A2'' E''", "E1/2 E3/2 E5/2"),
    ("D5", "A1 A2 E1 E2", "E1/2 E3/2 1E5/2 2E5/2"),
    ("D6h", "A1g A2g B1g B2g E1g E2g A1u A2u B1u B2u E1u E2u", "E1/2g E3/2g E5/2g E1/2u E3/2u E5/2u"),
    ("D4d", "A1 A2 B1 B2 E1 E2 E3", "E1/2 E3/2 E5/2 E7/2"),
]  # fmt: skip


@pytest.mark.parametrize(("symbol", "ordinary", "projective"), EXAMPLES)
def test_canonical_symbols_match_the_spec_examples(
    symbol: str, ordinary: str, projective: str
) -> None:
    assert symbols_of(symbol) == (ordinary, projective)


@pytest.mark.parametrize(
    ("symbol", "element", "expected"),
    [
        ("D4", "C2'(1)", {"B1": "1", "B2": "-1"}),
        ("C4v", "sigma_v(1)", {"B1": "1", "B2": "-1"}),
        ("D2d", "C2'(1)", {"B1": "1", "B2": "-1"}),
        ("D2", "C2", {"B1": "1", "B2": "-1", "B3": "-1"}),
        ("D2", "C2'(1)", {"B1": "-1", "B2": "-1", "B3": "1"}),
        ("D2", "C2''(1)", {"B1": "-1", "B2": "1", "B3": "-1"}),
        ("D3h", "S3^-1", {"E1/2": "sqrt3", "E3/2": "0", "E5/2": "-sqrt3"}),
        ("D2d", "S4^-1", {"E1/2": "sqrt2", "E3/2": "-sqrt2"}),
        ("D3", "C2'(1)", {"1E3/2": "i", "2E3/2": "-i"}),
        ("D3", "C2'(2)", {"1E3/2": "-i", "2E3/2": "i"}),
        ("C3", "C3^1", {"1E": "zeta", "2E": "zeta^2", "1E1/2": "-zeta^2", "2E1/2": "-zeta", "B3/2": "-1"}),
        ("D3h", "sigma_h", {"A1'": "1", "A1''": "-1", "E'": "2", "E''": "-2"}),
    ],
)  # fmt: skip
def test_characters_on_named_elements(
    symbol: str, element: str, expected: dict[str, str]
) -> None:
    document = build_document(Request("characters", symbol, values=("pretty",)))
    items = document["elements"]["items"]
    index = next(r["index"] for r in items if r["label"] == element)
    table = document["characters"]
    column = next(c["index"] for c in table["classes"] if index in c["members"])
    got = {
        r["symbol"]: table["pretty"][k][column] for k, r in enumerate(table["irreps"])
    }
    for symbol_text, value in expected.items():
        assert got[symbol_text] == value, (symbol_text, got)


def test_canonical_order_examples() -> None:
    assert symbols_of("D4h")[0].split()[:5] == ["A1g", "A2g", "B1g", "B2g", "Eg"]
    assert symbols_of("C6") == (
        "A B 1E1 2E1 1E2 2E2",
        "1E1/2 2E1/2 1E3/2 2E3/2 1E5/2 2E5/2",
    )
    assert symbols_of("D3h")[1] == "E1/2 E3/2 E5/2"


def test_symbol_sort_key_letter_before_subscript() -> None:
    keys = [
        Symbol.parse(text).sort_key(True)
        for text in ("B5/2", "1E1/2", "2E1/2", "1E3/2")
    ]
    assert keys == sorted(keys)
    assert Symbol.parse("1E3/2u").text == "1E3/2u"
    assert Symbol.parse("1E3/2u").partner == Symbol("E", "3/2", "u", 2)
    with pytest.raises(ValueError):
        Symbol.parse("F3/2")


@pytest.mark.parametrize("symbol", AXIAL_SEED)
def test_irrep_records_are_structured(symbol: str) -> None:
    table = characters_of(symbol)
    for record in table["irreps"]:
        parsed = Symbol.parse(record["symbol"])
        assert record["letter"] == parsed.letter
        assert record.get("subscript") == parsed.subscript
        assert record.get("suffix") == parsed.suffix
        assert record.get("component") == parsed.component
        if "partner" in record:
            partner = table["irreps"][record["partner"]]
            assert partner["partner"] == record["index"]
            assert partner["symbol"] == parsed.partner.text if parsed.partner else False
        assert record["frobenius_schur"] in (1, 0, -1)
        assert list(record) == [key for key in (
            "index", "symbol", "dimension", "projective", "letter", "component", "partner", "subscript", "suffix", "frobenius_schur",
        ) if key in record]  # fmt: skip


# ---- --irreps and --values ------------------------------------------------------------------


def test_irreps_filter_keeps_global_indices() -> None:
    projective = build_document(Request("characters", "C3v", irreps="projective"))
    assert [r["index"] for r in projective["characters"]["irreps"]] == [3, 4, 5]
    assert projective["irreps"] == "projective"
    assert projective["characters"]["class_group"] == "G_double"
    assert len(projective["characters"]["classes"]) == 6
    ordinary = build_document(Request("characters", "C3v", irreps="ordinary"))
    assert [r["index"] for r in ordinary["characters"]["irreps"]] == [0, 1, 2]
    assert ordinary["characters"]["class_group"] == "G"
    assert ordinary["characters"]["classes"] == ordinary["elements"]["classes"]
    everything = build_document(Request("characters", "C3v"))
    assert (
        everything["characters"]["irreps"]
        == ordinary["characters"]["irreps"] + projective["characters"]["irreps"]
    )
    assert (
        everything["characters"]["exact"]
        == everything["characters"]["exact"][:3] + projective["characters"]["exact"]
    )
    assert (
        everything["characters"]["classes"]
        == model_document("C3v", double=True)["elements"]["classes"]
    )
    with pytest.raises(UsageError):
        build_document(Request("characters", "C3v", irreps="spinor"))


def test_ordinary_rows_over_double_classes_repeat_on_the_barred_class() -> None:
    document = build_document(Request("characters", "D3"))
    table = document["characters"]
    columns = table["classes"]
    for k, record in enumerate(table["irreps"]):
        for c, cls in enumerate(columns):
            if cls["split"] and cls["barred"]:
                lifted = next(
                    c2
                    for c2, other in enumerate(columns)
                    if other["base_class"] == cls["base_class"] and not other["barred"]
                )
                sign = -1 if record["projective"] else 1
                assert table["numeric"][k][c] == [
                    sign * x or 0.0 for x in table["numeric"][k][lifted]
                ]


@pytest.mark.parametrize(
    "text",
    [
        "exact",
        "numeric",
        "pretty",
        "both",
        "exact,pretty",
        "numeric,pretty",
        "both,pretty",
        "pretty,exact,numeric",
    ],
)
def test_values_option_selects_the_parallel_arrays(text: str) -> None:
    forms = parse_values(text)
    document = build_document(Request("characters", "D3h", values=forms))
    assert document["values"] == list(forms)
    assert [form for form in VALUE_FORMS if form in forms] == list(forms)
    table = document["characters"]
    present = [form for form in VALUE_FORMS if form in table]
    assert present == list(forms)
    shapes = {form: [len(row) for row in table[form]] for form in present}
    assert len(set(map(tuple, shapes.values()))) == 1
    assert all(
        len(rows) == len(table["irreps"]) for rows in (table[form] for form in present)
    )


# ---- the prototype samples ----------------------------------------------------------------


@pytest.mark.parametrize(
    ("sample", "invocation"),
    [
        ("C3v_characters.json", Request("characters", "C3v")),
        (
            "C3v_characters_ordinary.json",
            Request("characters", "C3v", irreps="ordinary"),
        ),
        ("D6h_characters.json", Request("characters", "D6h")),
    ],
)
def test_characters_equal_the_prototype_samples(
    sample: str, invocation: Request
) -> None:
    path = SAMPLES / sample
    if not path.is_file():
        pytest.skip(f"the prototype sample {path.name} is not checked in")
    expected = json.loads(path.read_text())["characters"]
    table = build_document(invocation)["characters"]
    assert table["class_group"] == expected["class_group"]
    assert table["irreps"] == expected["irreps"]
    assert table["classes"] == expected["classes"]
    assert table["exact"] == expected["exact"]
    # the prototype rounded with plain complex arithmetic; the numeric forms agree to the noise
    for row, expected_row in zip(table["numeric"], expected["numeric"], strict=True):
        for value, expected_value in zip(row, expected_row, strict=True):
            assert value == pytest.approx(expected_value, abs=1e-15)


# ---- the schema (OS-1) and the polyhedral groups ---------------------------------------------


@pytest.mark.parametrize("which", ["ordinary", "projective", "all"])
@pytest.mark.parametrize("symbol", AXIAL_SEED)
def test_os1_characters_documents_validate(symbol: str, which: str) -> None:
    validator = jsonschema.Draft202012Validator(SCHEMA)
    for forms in (("exact", "numeric"), ("pretty",), ("exact", "numeric", "pretty")):
        document = build_document(
            Request("characters", symbol, values=forms, irreps=which)
        )
        validator.validate(document)
        assert document["tables"] == ["elements", "characters"]
        assert document["irreps"] == which


@pytest.mark.parametrize("symbol", POLYHEDRAL_GROUPS)
def test_polyhedral_characters_are_not_available_yet(symbol: str) -> None:
    assert "characters" not in model(symbol)
    with pytest.raises(PointGroupTablesError, match="not available"):
        build_document(Request("characters", symbol))


# ---- the field of the closed forms -----------------------------------------------------------


@pytest.mark.parametrize(
    "symbol", [s for s in DEFAULT_GROUPS if s not in POLYHEDRAL_GROUPS]
)
def test_canonical_conductor_is_minimal_for_the_closed_form_matrices(
    symbol: str,
) -> None:
    """``VL-3`` ahead of the matrices table: no maximal proper divisor of the canonical
    conductor holds every entry of every irrep matrix in the canonical basis."""
    computed = build_model(parse_group(symbol))
    F = computed.field
    M = F.conductor
    entries = list(
        itertools.chain.from_iterable(
            itertools.chain.from_iterable(matrix)
            for irrep in computed.irreps
            for matrix in irrep.matrices
        )
    )
    for p in range(2, M + 1):
        if M % p or any(p % q == 0 for q in range(2, p)):
            continue
        smaller = cyclotomic_field(M // p)
        assert any(F.descend(x, smaller) is None for x in entries), (symbol, M, p)
    assert M % 4 != 2


@pytest.mark.parametrize("symbol", AXIAL_SEED)
def test_every_character_lies_in_the_declared_field_and_traces_the_matrices(
    symbol: str,
) -> None:
    computed = build_model(parse_group(symbol))
    F = computed.field
    for irrep in computed.irreps:
        assert len(irrep.characters) == computed.order
        for matrix, character in zip(irrep.matrices, irrep.characters, strict=True):
            assert len(matrix) == irrep.dimension
            assert F.total(matrix[i][i] for i in range(irrep.dimension)) == character
