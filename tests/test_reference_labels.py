"""The canonical labels behind the Altmann & Herzig label transcriptions.

Every transcription under ``tests/reference/altmann-herzig/labels/`` carries a
``[map]`` from the book's labels to canonical labels; its values must be
exactly the canonical labels the tool emits for that group. This is the
element-identification half of ``REF-6``; the convention itself is a later
ticket.
"""

from __future__ import annotations

import tomllib
from pathlib import Path

import pytest
from conftest import model

LABELS = Path(__file__).resolve().parent / "reference" / "altmann-herzig" / "labels"
FILES = sorted(LABELS.glob("*.toml")) if LABELS.is_dir() else []


@pytest.mark.reference
@pytest.mark.parametrize("path", FILES, ids=lambda path: path.stem)
def test_transcription_map_uses_the_canonical_labels(path: Path) -> None:
    transcription = tomllib.loads(path.read_text())
    if "map" not in transcription:
        pytest.skip(f"{path.name} has no [map] table")
    ours = {
        record["label"] for record in model(transcription["group"])["elements"]["items"]
    }
    assert set(transcription["map"].values()) == ours


@pytest.mark.reference
def test_transcriptions_are_present() -> None:
    if not FILES:
        pytest.skip(f"no transcription under {LABELS}")
