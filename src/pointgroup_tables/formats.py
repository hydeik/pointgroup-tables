"""Serialisation of a document: the fixed JSON layout (docs/spec/output-schema.md,
"JSON formatting") and YAML and TOML as plain dumps of the same object
(docs/spec/values.md, "YAML and TOML").

The document model holds objects with string keys, arrays, integers, strings,
doubles and booleans, in a fixed key order and without ``None``. The three
emitters share one rule for arrays: an array of scalars (recursively, of
arrays of scalars) is written on one line when it fits in ``LINE_WIDTH``
columns and one entry per line otherwise. Every emitter writes ASCII only,
with LF line ends and a trailing newline.
"""

from __future__ import annotations

import json
import math
import re
from collections.abc import Callable
from typing import Any

LINE_WIDTH = 100


def _normalised(value: Any) -> Any:
    """Negative zero prints as ``0.0``; everything else is untouched."""
    if isinstance(value, list):
        return [_normalised(item) for item in value]
    if isinstance(value, float) and value == 0:
        return 0.0
    return value


def _is_leafy(value: Any) -> bool:
    """A scalar, or an array whose entries are (recursively) leafy."""
    if isinstance(value, dict):
        return False
    if isinstance(value, list):
        return all(_is_leafy(item) for item in value)
    return True


def _fits(text: str, prefix: int) -> bool:
    """Whether ``text`` fits on a line after ``prefix`` columns, with a trailing comma."""
    return prefix + len(text) + 1 <= LINE_WIDTH


# ---- JSON ---------------------------------------------------------------------


def _inline(value: Any) -> str:
    return json.dumps(
        _normalised(value), separators=(", ", ": "), ensure_ascii=True, allow_nan=False
    )


def _format(value: Any, indent: int, prefix: int) -> str:
    """Format ``value`` at ``indent`` columns; ``prefix`` columns precede it on its line."""
    pad = " " * indent
    if isinstance(value, dict):
        if not value:
            return "{}"
        lines = [
            f'{pad}  "{key}": {_format(item, indent + 2, indent + 6 + len(key))}'
            for key, item in value.items()
        ]
        return "{\n" + ",\n".join(lines) + f"\n{pad}}}"
    if isinstance(value, list):
        if not value:
            return "[]"
        if _is_leafy(value):
            text = _inline(value)
            if _fits(text, prefix):
                return text
        lines = [f"{pad}  {_format(item, indent + 2, indent + 2)}" for item in value]
        return "[\n" + ",\n".join(lines) + f"\n{pad}]"
    return _inline(value)


def to_json(document: dict[str, Any]) -> str:
    """The fixed JSON layout: one key per line, leaf arrays inline when they fit."""
    return _format(document, 0, 0) + "\n"


# ---- scalars shared by YAML and TOML -----------------------------------------

_SHORT_ESCAPES = {
    '"': '\\"',
    "\\": "\\\\",
    "\b": "\\b",
    "\f": "\\f",
    "\n": "\\n",
    "\r": "\\r",
    "\t": "\\t",
}


def _quoted(text: str) -> str:
    """A double-quoted string in the escape syntax YAML and TOML share, ASCII only."""
    parts = []
    for char in text:
        if char in _SHORT_ESCAPES:
            parts.append(_SHORT_ESCAPES[char])
        elif " " <= char <= "~":
            parts.append(char)
        elif ord(char) <= 0xFFFF:
            parts.append(f"\\u{ord(char):04X}")
        else:
            parts.append(f"\\U{ord(char):08X}")
    return '"' + "".join(parts) + '"'


def _float(value: float) -> str:
    """The shortest round-tripping form, spelt so that YAML 1.1 and TOML read a float:
    a dot in the mantissa and a signed exponent."""
    if not math.isfinite(value):
        raise ValueError(f"{value!r} is not a document value")
    text = repr(0.0 if value == 0 else value)
    if "e" in text:
        mantissa, exponent = text.split("e")
        if "." not in mantissa:
            mantissa += ".0"
        text = f"{mantissa}e{exponent}"
    return text


def _scalar(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return _float(value)
    if isinstance(value, str):
        return _quoted(value)
    raise TypeError(f"{type(value).__name__} is not a document value")


def _flow(value: Any) -> str:
    """A leafy value on one line: ``[a, b, [c]]`` reads the same in YAML and TOML."""
    if isinstance(value, list):
        return "[" + ", ".join(_flow(item) for item in value) + "]"
    return _scalar(value)


def _leaf_or_none(value: Any, prefix: int) -> str | None:
    """The one-line form of ``value`` when it is a scalar, empty, or a leafy array that fits."""
    if isinstance(value, dict):
        return "{}" if not value else None
    if isinstance(value, list):
        if not value:
            return "[]"
        if _is_leafy(value):
            text = _flow(value)
            return text if _fits(text, prefix) else None
        return None
    return _scalar(value)


# ---- YAML ---------------------------------------------------------------------


def _yaml_mapping(mapping: dict[str, Any], indent: int) -> list[str]:
    """The block lines of a non-empty mapping whose keys start at ``indent``."""
    pad = " " * indent
    lines: list[str] = []
    for key, item in mapping.items():
        head = f"{pad}{_quoted(key)}:"
        leaf = _leaf_or_none(item, len(head) + 1)
        if leaf is None:
            lines.append(head)
            lines.extend(_yaml_block(item, indent + 2))
        else:
            lines.append(f"{head} {leaf}")
    return lines


def _yaml_sequence(sequence: list[Any], indent: int) -> list[str]:
    """The block lines of a non-empty sequence whose dashes start at ``indent``."""
    pad = " " * indent
    lines: list[str] = []
    for item in sequence:
        leaf = _leaf_or_none(item, indent + 2)
        if leaf is not None:
            lines.append(f"{pad}- {leaf}")
        elif isinstance(item, dict):
            block = _yaml_mapping(item, indent + 2)
            lines.append(f"{pad}- {block[0][indent + 2 :]}")
            lines.extend(block[1:])
        else:
            lines.append(f"{pad}-")
            lines.extend(_yaml_sequence(item, indent + 2))
    return lines


def _yaml_block(value: dict[str, Any] | list[Any], indent: int) -> list[str]:
    if isinstance(value, dict):
        return _yaml_mapping(value, indent)
    return _yaml_sequence(value, indent)


def to_yaml(document: dict[str, Any]) -> str:
    """A block-style YAML dump: every key and string double-quoted (``n``, ``yes``
    and the like are booleans to a YAML 1.1 reader), leaf arrays in flow style."""
    if not document:
        return "{}\n"
    return "\n".join(_yaml_mapping(document, 0)) + "\n"


# ---- TOML ---------------------------------------------------------------------

_BARE_KEY = re.compile(r"[A-Za-z0-9_-]+")


def _toml_key(key: str) -> str:
    return key if _BARE_KEY.fullmatch(key) else _quoted(key)


def _toml_value(value: Any, indent: int, prefix: int) -> str:
    """A TOML value at ``indent`` columns; ``prefix`` columns precede it on its line."""
    pad = " " * indent
    if isinstance(value, dict):
        if not value:
            return "{}"
        pairs = []
        for key, item in value.items():
            name = _toml_key(key)
            pairs.append(
                f"{name} = {_toml_value(item, indent, prefix + len(name) + 3)}"
            )
        return "{ " + ", ".join(pairs) + " }"
    leaf = _leaf_or_none(value, prefix)
    if leaf is not None:
        return leaf
    lines = [f"{pad}  {_toml_value(item, indent + 2, indent + 2)}" for item in value]
    return "[\n" + ",\n".join(lines) + f"\n{pad}]"


def _toml_pairs(mapping: dict[str, Any]) -> list[str]:
    lines = []
    for key, item in mapping.items():
        name = _toml_key(key)
        lines.append(f"{name} = {_toml_value(item, 0, len(name) + 3)}")
    return lines


def to_toml(document: dict[str, Any]) -> str:
    """A TOML dump: the metadata keys first (scalars, arrays and inline tables),
    then one ``[table]`` section per table named in ``tables``; record lists are
    arrays of inline tables. Every other object is an inline table."""
    names = document.get("tables")
    sections = [
        key
        for key in (names if isinstance(names, list) else [])
        if isinstance(document.get(key), dict)
    ]
    lines: list[str] = []
    seen_section = False
    for key, item in document.items():
        if key in sections:
            seen_section = True
            lines.extend(["", f"[{_toml_key(key)}]", *_toml_pairs(item)])
        elif seen_section:
            raise ValueError(f"key '{key}' comes after a table section")
        else:
            lines.extend(_toml_pairs({key: item}))
    if lines and lines[0] == "":
        del lines[0]
    return "\n".join(lines) + "\n"


EMITTERS: dict[str, Callable[[dict[str, Any]], str]] = {
    "json": to_json,
    "yaml": to_yaml,
    "toml": to_toml,
}
