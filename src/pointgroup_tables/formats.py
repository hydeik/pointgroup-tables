"""Serialisation of a document (docs/spec/output-schema.md, "JSON formatting")."""

from __future__ import annotations

import json
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


def _inline(value: Any) -> str:
    return json.dumps(
        _normalised(value), separators=(", ", ": "), ensure_ascii=True, allow_nan=False
    )


def _is_leafy(value: Any) -> bool:
    """A scalar, or an array whose entries are (recursively) leafy."""
    if isinstance(value, dict):
        return False
    if isinstance(value, list):
        return all(_is_leafy(item) for item in value)
    return True


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
            if prefix + len(text) + 1 <= LINE_WIDTH:  # + 1 for a trailing comma
                return text
        lines = [f"{pad}  {_format(item, indent + 2, indent + 2)}" for item in value]
        return "[\n" + ",\n".join(lines) + f"\n{pad}]"
    return _inline(value)


def to_json(document: dict[str, Any]) -> str:
    """The fixed JSON layout: one key per line, leaf arrays inline when they fit."""
    return _format(document, 0, 0) + "\n"


EMITTERS: dict[str, Callable[[dict[str, Any]], str]] = {"json": to_json}
