"""The document layout across every subcommand (docs/spec/output-schema.md
"Verification", ``OS-2`` and ``OS-3``): ``tables`` and ``values`` list exactly
what is present and every form of a table has the same shape; every index is
in range for its index space, ``items[i].index = i``, ``classes[c].index = c``,
and the rows of ``characters`` and the entries of ``matrices`` follow the
canonical irrep order."""

from __future__ import annotations

from typing import Any

import pytest
from conftest import SEED_GROUPS

from pointgroup_tables.document import (
    TABLES_OF_SUBCOMMAND,
    VALUE_FORMS,
    Request,
    build_document,
)

SUBCOMMANDS = tuple(TABLES_OF_SUBCOMMAND)
METADATA_KEYS = (
    "schema_version",
    "generator",
    "group",
    "tables",
    "double",
    "irreps",
    "values",
    "field",
    "label_convention",
    "basis_convention",
)

# the option variants of one document, beyond the default request of every subcommand
VARIANTS = [
    Request("elements", "D4h", double=True),
    Request("cayley", "D4h", double=True),
    Request("all", "D4h", double=True),
    Request("characters", "D4h", irreps="ordinary"),
    Request("matrices", "D4h", irreps="projective"),
    Request("all", "D4h", irreps="projective", values=VALUE_FORMS),
    Request("all", "C3v", values=("pretty",)),
    Request("characters", "C4", irreps="ordinary", values=("numeric",)),
    Request("all", "Th", double=True, irreps="ordinary"),
]


def _requests() -> list[Any]:
    defaults = [
        pytest.param(Request(subcommand, symbol), id=f"{symbol}-{subcommand}")
        for symbol in SEED_GROUPS
        for subcommand in SUBCOMMANDS
    ]
    variants = [
        pytest.param(
            request,
            id=f"{request.group}-{request.subcommand}-{request.irreps}"
            f"{'-double' if request.double else ''}-{'+'.join(request.values)}",
        )
        for request in VARIANTS
    ]
    return defaults + variants


def _character_shape(table: dict[str, Any], form: str) -> list[int]:
    return [len(row) for row in table[form]]


def _matrix_shape(entry: dict[str, Any], form: str) -> list[tuple[int, list[int]]]:
    return [(len(matrix), [len(row) for row in matrix]) for matrix in entry[form]]


# ---- OS-2: tables, values, shapes --------------------------------------------------------


@pytest.mark.parametrize("request_", _requests())
def test_os2_tables_and_values_list_exactly_what_is_present(request_: Request) -> None:
    document = build_document(request_)
    keys = list(document)
    tables = document["tables"]
    assert tables == list(TABLES_OF_SUBCOMMAND[request_.subcommand])
    assert tables[0] == "elements"
    # the metadata comes first, in the fixed order; then the tables, in table order
    assert keys[: -len(tables)] == [k for k in METADATA_KEYS if k in document]
    assert keys[-len(tables) :] == tables
    assert ("irreps" in document) == bool({"characters", "matrices"} & set(tables))

    forms = document["values"]
    assert forms == [form for form in VALUE_FORMS if form in forms]  # fixed order
    assert forms == list(request_.values)
    for name in ("elements", "cayley", "factor"):
        if name in document:
            assert not set(VALUE_FORMS) & set(document[name])
    if "characters" in document:
        table = document["characters"]
        assert [form for form in VALUE_FORMS if form in table] == forms
        shapes = {tuple(_character_shape(table, form)) for form in forms}
        assert len(shapes) == 1
        (shape,) = shapes
        assert shape == (len(table["classes"]),) * len(table["irreps"])
    if "matrices" in document:
        for entry in document["matrices"]["irreps"]:
            assert [form for form in VALUE_FORMS if form in entry] == forms
            shapes = {str(_matrix_shape(entry, form)) for form in forms}
            assert len(shapes) == 1


# ---- OS-3: index spaces and canonical order ---------------------------------------------


def _check_class_records(
    classes: list[dict[str, Any]], element_count: int, double: bool
) -> None:
    for c, record in enumerate(classes):
        assert record["index"] == c
        members = record["members"]
        assert 0 <= record["representative"] < element_count
        assert (
            members[0] == record["representative"]
        )  # representative first, rest ascending
        assert members[1:] == sorted(members[1:])
        assert all(0 <= member < element_count for member in members)
        assert len(set(members)) == len(members) == record["size"]
        assert 0 <= record["inverse"] < len(classes)
        assert classes[record["inverse"]]["inverse"] == c
        if double:
            assert "regular" not in record
            assert 0 <= record["base_class"] <= c
            assert isinstance(record["split"], bool)
            # ``barred`` is carried by split classes only
            assert ("barred" in record) == record["split"]
            assert record.get("barred", False) in (True, False)
        else:
            assert isinstance(record["regular"], bool)
            assert not {"base_class", "split", "barred"} & set(record)


@pytest.mark.parametrize("request_", _requests())
def test_os3_every_index_is_in_range_for_its_index_space(request_: Request) -> None:
    document = build_document(request_)
    order = document["group"]["order"]
    elements = document["elements"]
    items = elements["items"]
    assert len(items) == (2 * order if document["double"] else order)
    for i, item in enumerate(items):
        assert item["index"] == i
        assert 0 <= item["inverse"] < len(items)
        assert items[item["inverse"]]["inverse"] == i
        assert 0 <= item["class"] < len(elements["classes"])
        assert i in elements["classes"][item["class"]]["members"]
    _check_class_records(elements["classes"], len(items), document["double"])
    if "cayley" in document:
        table = document["cayley"]["table"]
        assert len(table) == len(items)
        assert all(len(row) == len(items) for row in table)
        assert all(0 <= entry < len(items) for row in table for entry in row)
    if "factor" in document:
        table = document["factor"]["table"]  # always on G
        assert len(table) == order
        assert all(len(row) == order for row in table)
        assert all(entry in (1, -1) for row in table for entry in row)
    if "characters" in document:
        table = document["characters"]
        indices = [record["index"] for record in table["irreps"]]
        assert indices == sorted(set(indices))
        assert [record["projective"] for record in table["irreps"]] == sorted(
            record["projective"] for record in table["irreps"]
        )
        # the barred rule: an index i >= |G| names the barred element of i - |G|
        space = order if table["class_group"] == "G" else 2 * order
        _check_class_records(table["classes"], space, table["class_group"] != "G")
    if "matrices" in document:
        entries = document["matrices"]["irreps"]
        indices = [entry["index"] for entry in entries]
        assert indices == sorted(set(indices))
        for entry in entries:
            for form in document["values"]:
                assert len(entry[form]) == order  # always indexed by the elements of G
    if "characters" in document and "matrices" in document:
        assert [e["index"] for e in document["matrices"]["irreps"]] == [
            r["index"] for r in document["characters"]["irreps"]
        ]


def test_os3_irrep_indices_are_global_and_in_canonical_order() -> None:
    """Filtering removes rows but does not renumber; the rows of a filtered table are
    the rows of the unfiltered one with the same records."""
    full = build_document(Request("characters", "C3v"))["characters"]
    assert [r["index"] for r in full["irreps"]] == list(range(len(full["irreps"])))
    for which in ("ordinary", "projective"):
        table = build_document(Request("characters", "C3v", irreps=which))["characters"]
        kept = [r for r in full["irreps"] if r["projective"] == (which == "projective")]
        assert table["irreps"] == kept
