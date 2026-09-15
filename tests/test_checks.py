"""The checkers detect corrupted models (a first slice of the mutation tests)."""

from __future__ import annotations

import copy
from typing import Any

from conftest import model

from pointgroup_tables import checks


def corrupted(symbol: str) -> dict[str, Any]:
    return copy.deepcopy(model(symbol))


def failed(document: dict[str, Any]) -> set[str]:
    return {failure.check for failure in checks.check(document)}


def test_ag1_catches_a_replaced_cayley_entry() -> None:
    document = corrupted("C3v")
    document["cayley"]["table"][1][2] = 3
    assert "AG-1" in failed(document)


def test_ag2_catches_a_missing_class() -> None:
    document = corrupted("Oh")
    document["elements"]["classes"].pop()
    assert "AG-2" in failed(document)


def test_cl1_catches_an_element_moved_to_another_class() -> None:
    document = corrupted("C3v")
    document["elements"]["classes"][1]["members"].append(3)
    document["elements"]["classes"][2]["members"].remove(3)
    assert "CL-1" in failed(document)


def test_cl2_catches_a_bad_or_duplicate_label() -> None:
    document = corrupted("C3v")
    document["elements"]["classes"][1]["label"] = "2C3x"
    assert failed(document) == {"CL-2"}
    document["elements"]["classes"][1]["label"] = "E"
    assert failed(document) == {"CL-2"}


def test_cl3_catches_a_flipped_regular_flag() -> None:
    for symbol, index in (("Oh", 1), ("C3v", 2), ("D3h", 3)):
        document = corrupted(symbol)
        record = document["elements"]["classes"][index]
        record["regular"] = not record["regular"]
        assert failed(document) == {"CL-3"}, symbol


def test_cl4_catches_a_wrong_inverse_class() -> None:
    document = corrupted("C6")
    document["elements"]["classes"][1]["inverse"] = 1
    assert failed(document) == {"CL-4"}


def test_cl5_catches_swapped_polyhedral_labels() -> None:
    document = corrupted("Td")
    classes = document["elements"]["classes"]
    classes[1]["label"], classes[2]["label"] = classes[2]["label"], classes[1]["label"]
    assert "CL-5" in failed(document)


def test_ag3_catches_a_reordered_parity_block() -> None:
    document = corrupted("C3v")
    document["elements"]["items"][0]["parity"] = -1
    assert failed(document) >= {"AG-3"}


def test_clean_models_pass_every_applicable_check() -> None:
    for symbol in ("C1", "C3v", "D6h", "Oh", "Ih"):
        assert checks.check(model(symbol)) == []
