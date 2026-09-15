"""The fixed JSON layout (docs/spec/output-schema.md, "JSON formatting")."""

from __future__ import annotations

import json

from pointgroup_tables.formats import to_json


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
