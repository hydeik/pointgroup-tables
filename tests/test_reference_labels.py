"""``REF-6``: the ``altmann-herzig`` label convention reproduces the transcribed
element labels of every group, element by element (docs/spec/label-conventions.md,
docs/spec/verification.md).

Every transcription under ``tests/reference/altmann-herzig/labels/`` carries a
``[map]`` from the book's labels to canonical labels. Its values must be
exactly the canonical labels the tool emits for that group (the element
identification), and the reference label the convention attaches to each
element must be the book's label of it.
"""

from __future__ import annotations

import tomllib
from pathlib import Path

import pytest
from conftest import model
from test_label_convention import AH

from pointgroup_tables.document import model_document

LABELS = Path(__file__).resolve().parent / "reference" / "altmann-herzig" / "labels"
FILES = sorted(LABELS.glob("*.toml")) if LABELS.is_dir() else []


def _transcription(path: Path) -> dict:
    transcription = tomllib.loads(path.read_text())
    if "map" not in transcription:
        pytest.skip(f"{path.name} has no [map] table")
    return transcription


@pytest.mark.reference
@pytest.mark.parametrize("path", FILES, ids=lambda path: path.stem)
def test_transcription_map_uses_the_canonical_labels(path: Path) -> None:
    transcription = _transcription(path)
    ours = {
        record["label"] for record in model(transcription["group"])["elements"]["items"]
    }
    assert set(transcription["map"].values()) == ours


@pytest.mark.reference
@pytest.mark.parametrize("path", FILES, ids=lambda path: path.stem)
def test_ref6_reference_labels_reproduce_the_book(path: Path) -> None:
    transcription = _transcription(path)
    book_of = {canonical: book for book, canonical in transcription["map"].items()}
    document = model_document(transcription["group"], label_convention=AH)
    assert document["label_convention"]["coverage"] == "reference"
    mismatches = {
        record["label"]: (record["reference_label"], book_of[record["label"]])
        for record in document["elements"]["items"]
        if record["reference_label"] != book_of[record["label"]]
    }
    assert mismatches == {}, "canonical: (ours, book)"


@pytest.mark.reference
def test_transcriptions_are_present() -> None:
    if not FILES:
        pytest.skip(f"no transcription under {LABELS}")
