"""The checkers detect corrupted models (a first slice of the mutation tests)."""

from __future__ import annotations

import copy
from typing import Any

from conftest import double_model, model

from pointgroup_tables import checks
from pointgroup_tables.field import cyclotomic_field


def corrupted(symbol: str) -> dict[str, Any]:
    return copy.deepcopy(model(symbol))


def applicable(document: dict[str, Any]) -> list[str]:
    """Every applicable check but the all-pairs ``MX-2`` on the icosahedral groups,
    which the ``slow`` battery runs (``tests/test_matrices.py``)."""
    names = checks.applicable_checks(document)
    if document["group"]["order"] > 48:
        names.remove("MX-2")
    return names


def failed(document: dict[str, Any]) -> set[str]:
    return {failure.check for failure in checks.check(document, applicable(document))}


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
        document = model(symbol)
        assert checks.check(document, applicable(document)) == []


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
        document = double_model(symbol)
        assert checks.check(document, applicable(document)) == []


# ---- IL-*, CH-* and VL-* --------------------------------------------------------


def test_mut4_a_replaced_character_is_caught_by_ch1_and_ch2() -> None:
    """``MUT-4``: one character replaced by another value of the field."""
    for symbol in ("C3v", "D6h"):
        document = corrupted(symbol)
        table = document["characters"]
        row, column = next(
            (k, c)
            for k, values in enumerate(table["exact"])
            for c, value in enumerate(values)
            if value != values[0]
        )
        table["exact"][row][column] = table["exact"][row][
            0
        ]  # the dimension, another value of the field
        failures = failed(document)
        assert {"CH-1", "CH-2"} <= failures, symbol


def test_mut5_swapped_irrep_symbols_are_caught_by_il2() -> None:
    """``MUT-5``: two irrep symbols swapped."""
    document = corrupted("C3v")
    irreps = document["characters"]["irreps"]
    for key in ("symbol", "letter", "subscript"):
        irreps[0][key], irreps[1][key] = irreps[1].get(key), irreps[0].get(key)
    assert "IL-2" in failed(document)
    document = corrupted("D3")
    irreps = document["characters"]["irreps"]
    irreps[4]["symbol"], irreps[5]["symbol"] = irreps[5]["symbol"], irreps[4]["symbol"]
    irreps[4]["component"], irreps[5]["component"] = 2, 1
    irreps[4]["partner"], irreps[5]["partner"] = 5, 4
    assert {"IL-2", "IL-4"} <= failed(document)
    # the polyhedral rules: T1 / T2 of O by the character on C4^1(1), E1/2 / E7/2
    # of I by the character on s(C5^1(1))
    for symbol, first, second in (
        ("O", "T1", "T2"),
        ("I", "E1/2", "E7/2"),
        ("Th", "1Eg", "2Eg"),
    ):
        document = corrupted(symbol)
        irreps = document["characters"]["irreps"]
        a = next(r for r in irreps if r["symbol"] == first)
        b = next(r for r in irreps if r["symbol"] == second)
        for key in ("symbol", "subscript", "component"):
            a[key], b[key] = b.get(key), a.get(key)
            for r in (a, b):
                if r[key] is None:
                    del r[key]
        assert "IL-2" in failed(document), symbol


def test_il1_catches_a_bad_duplicate_or_inconsistent_symbol() -> None:
    document = corrupted("C3v")
    document["characters"]["irreps"][0]["symbol"] = "F1"
    assert "IL-1" in failed(document)
    document = corrupted("C3v")
    document["characters"]["irreps"][0]["symbol"] = "A2"
    assert "IL-1" in failed(document)
    document = corrupted("C3v")
    document["characters"]["irreps"][0]["letter"] = "B"
    assert "IL-1" in failed(document)


def test_il3_catches_a_wrong_dimension_indicator_or_partner() -> None:
    document = corrupted("C3v")
    document["characters"]["irreps"][2]["dimension"] = 1
    assert "IL-3" in failed(document)
    document = corrupted("C3v")
    document["characters"]["irreps"][3]["frobenius_schur"] = 1
    assert {"IL-3", "CH-5"} <= failed(document)
    document = corrupted("C3v")
    document["characters"]["irreps"][4]["partner"] = 4
    assert "IL-3" in failed(document)
    document = corrupted("C3v")
    del document["characters"]["irreps"][4]["partner"]
    del document["characters"]["irreps"][4]["component"]
    assert "IL-3" in failed(document)


def test_ch3_catches_a_missing_row() -> None:
    document = corrupted("C3v")
    table = document["characters"]
    for key in ("irreps", "exact", "numeric", "pretty"):
        table[key].pop()
    assert "CH-3" in failed(document)


def test_ch4_and_ch6_catch_a_projective_value_on_a_merged_class() -> None:
    document = corrupted("D4")
    table = document["characters"]
    merged = next(c for c, cls in enumerate(table["classes"]) if not cls["split"])
    projective = next(k for k, r in enumerate(table["irreps"]) if r["projective"])
    table["exact"][projective][merged] = table["exact"][0][0]
    assert {"CH-4", "CH-6"} <= failed(document)


def test_ch4_catches_a_barred_class_that_does_not_repeat_the_lifted_one() -> None:
    document = corrupted("C3v")
    table = document["characters"]
    barred = next(
        c for c, cls in enumerate(table["classes"]) if cls["split"] and cls["barred"]
    )
    table["exact"][0][barred] = table["exact"][1][
        barred
    ]  # A1 stays 1 on ~E: swap with A2's -> still 1
    table["exact"][2][barred] = table["exact"][0][
        0
    ]  # E on the barred class of E: 1 instead of 2
    assert "CH-4" in failed(document)


def test_vl1_and_vl2_catch_bad_value_forms() -> None:
    document = corrupted("C3v")
    document["characters"]["exact"][0][0] = [[2, 2], [0, 1], [0, 1], [0, 1]]
    assert "VL-1" in failed(document)
    document = corrupted("C3v")
    document["characters"]["exact"][0][0] = [[1, 1], [0, 1]]
    assert "VL-1" in failed(document)
    document = corrupted("C3v")
    document["characters"]["numeric"][0][0] = [1.0000000000000002, 0.0]
    assert "VL-2" in failed(document)
    document = corrupted("C3v")
    document["characters"]["numeric"][0][0] = [1.0, -0.0]
    assert "VL-2" in failed(document)
    document = corrupted("C3v")
    document["characters"]["pretty"][0][0] = "2"
    assert "VL-2" in failed(document)
    document = corrupted("C3v")
    document["characters"]["pretty"][0][0] = "sqrt2"
    assert "VL-2" in failed(document)


def test_checks_report_a_characters_table_without_an_exact_form() -> None:
    document = corrupted("C3v")
    del document["characters"]["exact"]
    assert {"IL-2", "CH-1"} <= failed(document)


# ---- the matrices table: MUT-2, MUT-3 and the other matrix checkers -----------------


def test_mut2_a_flipped_factor_sign_is_caught_by_dg2_and_mx2() -> None:
    document = corrupted("C3v")
    document["factor"]["table"][1][2] *= -1
    assert {"DG-2", "MX-2"} <= failed(document)


def test_mut3_a_matrix_multiplied_by_a_phase_is_caught_by_mx1_or_mx2_and_mx4() -> None:
    document = corrupted("C3v")
    F = cyclotomic_field(document["field"]["conductor"])
    entry = document["matrices"]["irreps"][2]  # E
    entry["exact"][1] = [
        [F.exact(F.mul(F.zeta(1), F.from_exact(pairs))) for pairs in row]
        for row in entry["exact"][1]
    ]
    names = failed(document)
    assert names & {"MX-1", "MX-2"} and "MX-4" in names


def test_mx3_catches_a_matrix_that_is_not_unitary() -> None:
    document = corrupted("C3v")
    entry = document["matrices"]["irreps"][2]
    entry["exact"][1] = [
        [[[2 * p, q] for p, q in pairs] for pairs in row] for row in entry["exact"][1]
    ]
    assert "MX-3" in failed(document)


def test_mx5_catches_a_reordered_basis() -> None:
    # swapping the two basis vectors of E keeps a unitary homomorphism with the
    # same traces, and only the eigenvalue order of D(c) tells
    document = corrupted("C3v")
    entry = document["matrices"]["irreps"][2]
    for form in ("exact", "numeric", "pretty"):
        entry[form] = [
            [[matrix[1][1], matrix[1][0]], [matrix[0][1], matrix[0][0]]]
            for matrix in entry[form]
        ]
    assert failed(document) == {"MX-5"}


def test_mx7_catches_a_missing_matrix_and_a_foreign_entry() -> None:
    document = corrupted("C3v")
    document["matrices"]["irreps"][0]["exact"].pop()
    assert "MX-7" in failed(document)
    document = corrupted("C3v")
    document["matrices"]["irreps"][0]["exact"][0][0][0] = [
        [1, 1],
        [0, 1],
    ]  # degree 2, not 4
    assert "MX-7" in failed(document)


def test_mx8_catches_a_duplicated_irrep() -> None:
    document = corrupted("C3v")
    entries = document["matrices"]["irreps"]
    entries[1]["exact"] = copy.deepcopy(entries[0]["exact"])
    entries[1]["numeric"] = copy.deepcopy(entries[0]["numeric"])
    assert {"MX-4", "MX-8"} <= failed(document)
