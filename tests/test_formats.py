"""The emitters: the fixed JSON layout, and YAML and TOML as plain dumps of the
same object (docs/spec/output-schema.md, "JSON formatting"; docs/spec/values.md,
"YAML and TOML")."""

from __future__ import annotations

import json
import tomllib

import pytest
import yaml
from conftest import SEED_GROUPS

from pointgroup_tables.document import Request, build_document
from pointgroup_tables.formats import EMITTERS, to_json, to_toml, to_yaml


def test_objects_one_key_per_line_and_leaf_arrays_inline() -> None:
    document = {
        "a": 1,
        "b": [1, 2, [3, 4]],
        "c": {"d": [], "e": {}},
        "f": [{"g": 1}],
        "h": "xé",
    }
    assert to_json(document) == (
        "{\n"
        '  "a": 1,\n'
        '  "b": [1, 2, [3, 4]],\n'
        '  "c": {\n'
        '    "d": [],\n'
        '    "e": {}\n'
        "  },\n"
        '  "f": [\n'
        "    {\n"
        '      "g": 1\n'
        "    }\n"
        "  ],\n"
        '  "h": "x\\u00e9"\n'
        "}\n"
    )


def test_long_leaf_arrays_break_one_entry_per_line() -> None:
    document = {"row": list(range(60))}
    text = to_json(document)
    assert text.startswith('{\n  "row": [\n    0,\n    1,\n')
    assert json.loads(text) == document


def test_floats_and_negative_zero() -> None:
    assert (
        to_json({"x": [0.5, -0.0, 1e-7, 2.0]})
        == '{\n  "x": [0.5, 0.0, 1e-07, 2.0]\n}\n'
    )


# ---- YAML and TOML: plain dumps of the same object --------------------------

TRICKY = {
    "n": 3,  # `n` unquoted is a YAML 1.1 boolean
    "yes": "no",  # so are `yes` and `no`
    "null": "~",
    "number_like": "1e3",
    "text": 'a "quoted" \\ back\\slash\n\ttab xé \U0001f600',
    "flag": True,
    "off": False,
    "floats": [0.5, -0.0, 1e-07, 2.0, 1e16, -3.25e-300, 123456789.0],
    "empty_list": [],
    "empty_dict": {},
    "nested": {"deep": {"deeper": [[1, 2], [3, [4, 5]]]}},
    "records": [
        {"index": 0, "members": [0, 1], "inner": {"k": "v"}},
        {"index": 1, "members": [], "inner": {}},
    ],
    "long": list(range(60)),
    "matrix": [[list(range(30)), list(range(30))], [[1], [2]]],
    "list_of_lists_of_records": [[{"a": 1}], [{"b": 2}, {"c": [{"d": 3}]}]],
    "list_of_lists": [[[[1]]], [[[]]]],
    "wide_record": [{"m": [list(range(40)), list(range(40))], "z": 0}],
}


def _expected(value: object) -> object:
    """The document as a parser returns it: negative zero is plain zero."""
    if isinstance(value, dict):
        return {key: _expected(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_expected(item) for item in value]
    if isinstance(value, float) and value == 0:
        return 0.0
    return value


def test_yaml_round_trips_through_pyyaml() -> None:
    text = to_yaml(TRICKY)
    assert yaml.safe_load(text) == _expected(TRICKY)
    assert text.isascii() and text.endswith("\n") and "\r" not in text


def test_toml_round_trips_through_tomllib() -> None:
    text = to_toml(TRICKY)
    assert tomllib.loads(text) == _expected(TRICKY)
    assert text.isascii() and text.endswith("\n") and "\r" not in text


def test_yaml_layout() -> None:
    document = {"a": 1, "n": "x", "b": [1, [2, 3]], "c": {"d": []}, "e": [{"f": 1, "g": {"h": 2}}, {"i": [1]}], "j": [[{"k": 1}]]}  # fmt: skip
    assert to_yaml(document) == (
        '"a": 1\n'
        '"n": "x"\n'
        '"b": [1, [2, 3]]\n'
        '"c":\n'
        '  "d": []\n'
        '"e":\n'
        '  - "f": 1\n'
        '    "g":\n'
        '      "h": 2\n'
        '  - "i": [1]\n'
        '"j":\n'
        "  -\n"
        '    - "k": 1\n'
    )


def test_toml_layout_metadata_first_then_one_section_per_table() -> None:
    document = {
        "schema_version": 1,
        "generator": {"name": "t", "version": "0"},
        "tables": ["elements", "cayley"],
        "values": [],
        "field": {"x": 0.5, "y": {}},
        "elements": {"group": "G", "items": [{"index": 0, "members": [0]}, {"index": 1, "members": []}]},
        "cayley": {"table": [[0, 1], [1, 0]]},
    }  # fmt: skip
    assert to_toml(document) == (
        "schema_version = 1\n"
        'generator = { name = "t", version = "0" }\n'
        'tables = ["elements", "cayley"]\n'
        "values = []\n"
        "field = { x = 0.5, y = {} }\n"
        "\n"
        "[elements]\n"
        'group = "G"\n'
        "items = [\n"
        "  { index = 0, members = [0] },\n"
        "  { index = 1, members = [] }\n"
        "]\n"
        "\n"
        "[cayley]\n"
        "table = [[0, 1], [1, 0]]\n"
    )


def test_toml_long_arrays_break_one_entry_per_line() -> None:
    document = {"row": list(range(60)), "grid": [list(range(20))] * 3}
    text = to_toml(document)
    assert text.startswith("row = [\n  0,\n  1,\n")
    assert "\ngrid = [\n  [0, 1, " in text and text.endswith(
        "19],\n  [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]\n]\n"
    )
    assert tomllib.loads(text) == document


def test_toml_keys_that_are_not_bare_are_quoted() -> None:
    text = to_toml({"a b": 1, "": 2, "x.y": {"k": 3}, "s": [{"": 4}]})
    assert tomllib.loads(text) == {"a b": 1, "": 2, "x.y": {"k": 3}, "s": [{"": 4}]}


def test_toml_rejects_a_document_with_metadata_after_a_table() -> None:
    with pytest.raises(ValueError, match="after"):
        to_toml({"tables": ["t"], "t": {"a": 1}, "late": 2})


def test_emitters_are_registered_under_the_format_names() -> None:
    assert list(EMITTERS) == ["json", "yaml", "toml"]


@pytest.mark.parametrize("symbol", SEED_GROUPS)
def test_os_elements_round_trips_through_yaml_and_toml(symbol: str) -> None:
    document = build_document(
        Request("characters", symbol, values=("exact", "numeric", "pretty"))
    )
    reference = json.loads(to_json(document))
    assert reference == document
    assert yaml.safe_load(to_yaml(document)) == reference
    assert tomllib.loads(to_toml(document)) == reference
