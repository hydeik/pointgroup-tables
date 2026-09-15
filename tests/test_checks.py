"""The checkers detect corrupted models (a first slice of the mutation tests)."""

from __future__ import annotations

import copy
from typing import Any

from conftest import double_model, model

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
        assert failed(document) == {"CL-3", "DG-6"}, (
            symbol
        )  # DG-6 counts the regular classes


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


# ---- DG-* and the G-dagger side of CL-* -------------------------------------


def corrupted_double(symbol: str) -> dict[str, Any]:
    return copy.deepcopy(double_model(symbol))


def test_dg2_catches_a_flipped_factor_sign() -> None:
    """``MUT-2``: one ``omega(g, h)`` sign flipped is caught by ``DG-2``."""
    for symbol in ("C3v", "Oh"):
        document = corrupted(symbol)
        document["factor"]["table"][1][2] *= -1
        assert "DG-2" in failed(document), symbol
        document = corrupted_double(symbol)
        document["factor"]["table"][1][2] *= -1
        assert "DG-2" in failed(document), symbol


def test_dg1_catches_a_bad_value_and_a_denormalised_table() -> None:
    document = corrupted("C3v")
    document["factor"]["table"][0][1] = -1
    assert "DG-1" in failed(document)
    document = corrupted("C3v")
    document["factor"]["table"][2][3] = 0
    assert "DG-1" in failed(document)


def test_dg3_catches_an_asymmetric_inverse_pair() -> None:
    document = corrupted("C6")
    labels = [record["label"] for record in document["elements"]["items"]]
    c6, c6_inverse = labels.index("C6^1"), labels.index("C6^-1")
    document["factor"]["table"][c6][c6_inverse] *= -1
    assert "DG-3" in failed(document)


def test_dg4_catches_a_corrupted_double_group_table() -> None:
    document = corrupted_double("C3v")
    N = 6
    document["cayley"]["table"][N][1] = 1  # ~E no longer bars C3
    assert "DG-4" in failed(document)
    document = corrupted_double("C3v")
    document["elements"]["items"][7]["label"] = "~C3^-1"
    assert "DG-4" in failed(document)
    document = corrupted_double("C3v")
    document["elements"]["items"][7]["base_index"] = 2
    assert "DG-4" in failed(document)


def test_dg5_catches_a_binary_rotation_that_does_not_square_to_the_full_turn() -> None:
    document = corrupted("D2h")
    labels = [record["label"] for record in document["elements"]["items"]]
    c2 = labels.index("C2")
    document["factor"]["table"][c2][c2] = 1
    assert "DG-5" in failed(document)
    document = corrupted("D2h")
    sigma, c2 = labels.index("sigma_h"), labels.index("C2")
    row = document["factor"]["table"][sigma]
    row[c2] *= -1  # omega(sigma_h, C2) no longer equals omega(C2, C2)
    assert "DG-5" in failed(document)


def test_dg6_catches_a_missing_double_class_and_a_wrong_regular_count() -> None:
    document = corrupted_double("D3h")
    document["elements"]["classes"].pop()
    assert "DG-6" in failed(document)
    document = corrupted("D3h")
    document["elements"]["classes"][2]["regular"] = True
    assert "DG-6" in failed(document)


def test_cl3_catches_a_flipped_split_flag_of_the_double_group() -> None:
    document = corrupted_double("D4")
    classes = document["elements"]["classes"]
    merged = next(cls for cls in classes if not cls["split"])
    merged["split"] = True
    merged["barred"] = False
    assert "CL-3" in failed(document)
    document = corrupted_double("D4")
    split = next(cls for cls in document["elements"]["classes"] if cls["split"])
    split["barred"] = not split["barred"]
    assert "CL-3" in failed(document)


def test_cl4_catches_a_wrong_base_class_of_the_double_group() -> None:
    document = corrupted_double("C3v")
    classes = document["elements"]["classes"]
    classes[2]["base_class"], classes[3]["base_class"] = 0, 0
    assert "CL-4" in failed(document)
    document = corrupted_double("C3v")
    classes = document["elements"]["classes"]
    del classes[3]["barred"]
    assert "CL-4" in failed(document)


def test_clean_double_models_pass_every_applicable_check() -> None:
    for symbol in ("C1", "C3v", "D6h", "Oh", "Ih"):
        assert checks.check(double_model(symbol)) == []
