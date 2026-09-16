"""Helpers shared by the Altmann & Herzig reference tests and the convention tool
(``tools/ah_convention.py``): reading the character transcriptions and matching
their columns to the classes of a labelled document."""

from __future__ import annotations

from typing import Any

from pointgroup_tables.field import CyclotomicField, Value
from pointgroup_tables.values import parse

# The character transcription spells these class labels differently from the labels
# transcription (whose README fixes ``sigma_x``) and from the reference class label
# of the convention (``15sigma`` of ``Ih``, docs/spec/conjugacy-classes.md).
CLASS_SPELLING = {
    "sigmax": "sigma_x",
    "sigmay": "sigma_y",
    "sigmaz": "sigma_z",
    "15sigma_h": "15sigma",
}

BookRow = tuple[str, list[Value], bool]  # label, characters, projective


def book_rows(
    data: dict[str, Any], F_book: CyclotomicField, F: CyclotomicField
) -> list[BookRow]:
    """The rows of a character transcription, lifted to the field ``F``."""
    return [
        (
            row["label"],
            [F_book.lift(parse(F_book, text), F) for text in row["characters"]],
            row["projective"],
        )
        for row in data["irreps"]
    ]


def swap_projective(
    rows: list[BookRow], classes: list[str], a: str, b: str
) -> list[BookRow]:
    """The rows with the columns ``a`` and ``b`` of the projective rows interchanged
    (a ``[[known_discrepancy]]`` entry with ``swap_projective``)."""
    ia, ib = classes.index(a), classes.index(b)
    out = []
    for label, values, projective in rows:
        values = list(values)
        if projective:
            values[ia], values[ib] = values[ib], values[ia]
        out.append((label, values, projective))
    return out


def corrected_rows(
    data: dict[str, Any], F_book: CyclotomicField, F: CyclotomicField
) -> list[BookRow]:
    """The rows with every documented discrepancy applied."""
    rows = book_rows(data, F_book, F)
    labels = [c["label"] for c in data["classes"]]
    for entry in data.get("known_discrepancy", []):
        rows = swap_projective(rows, labels, *entry["swap_projective"])
    return rows


def column_of(columns: list[dict[str, Any]], book_label: str) -> int:
    """The class of ``G-dagger`` the book's column names: the class whose reference
    label is the book's (the untilded half of a split class) or the merged class
    whose unbarred half it is."""
    label = CLASS_SPELLING.get(book_label, book_label)
    matches = [
        c
        for c, column in enumerate(columns)
        if column["reference_label"] == label
        or column["reference_label"].startswith(label + "+")
    ]
    if len(matches) != 1:
        raise LookupError((book_label, [c["reference_label"] for c in columns]))
    return matches[0]
