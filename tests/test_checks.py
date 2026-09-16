"""The mutation tests (docs/spec/verification.md, ``MUT-1`` ... ``MUT-7``): a
corrupted document model is caught by the named checks, and no other checker
crashes on it. Every check ``--verify`` runs is the target of at least one
mutation here (the closing test enforces it)."""

from __future__ import annotations

import copy
from typing import Any

import pytest
from conftest import double_model, model

from pointgroup_tables import checks
from pointgroup_tables.field import cyclotomic_field

MUTATION_GROUPS = (
    "C3v",
    "Oh",
)  # one axial, one polyhedral, both with projective irreps

# the checks the mutations of this module have been asserted to catch
TARGETED: set[str] = set()


def corrupted(symbol: str) -> dict[str, Any]:
    return copy.deepcopy(model(symbol))


def corrupted_double(symbol: str) -> dict[str, Any]:
    return copy.deepcopy(double_model(symbol))


def applicable(document: dict[str, Any]) -> list[str]:
    """Every applicable check but the all-pairs ``MX-2`` on the icosahedral groups,
    which the ``slow`` battery runs (``tests/test_matrices.py``)."""
    names = checks.applicable_checks(document)
    if document["group"]["order"] > 48:
        names.remove("MX-2")
    return names


def failed(document: dict[str, Any]) -> set[str]:
    """The checks that fail on a document; every applicable checker runs, so a
    checker that crashes on the corruption fails the test."""
    return {failure.check for failure in checks.check(document, applicable(document))}


def expect(
    document: dict[str, Any],
    *names: str,
    any_of: tuple[str, ...] = (),
    only: bool = False,
    note: object = None,
) -> set[str]:
    """Assert that the mutation is caught by every check in ``names`` (and by at
    least one of ``any_of``; by nothing else with ``only``), and record the targets."""
    names_failed = failed(document)
    assert set(names) <= names_failed, (names_failed, note)
    if any_of:
        assert names_failed & set(any_of), (names_failed, note)
    if only:
        assert names_failed == set(names), (names_failed, note)
    TARGETED.update(names)
    TARGETED.update(names_failed & set(any_of))
    return names_failed


def test_ag1_catches_a_replaced_cayley_entry() -> None:
    document = corrupted("C3v")
    document["cayley"]["table"][1][2] = 3
    expect(document, "AG-1")


def test_ag2_catches_a_missing_class() -> None:
    document = corrupted("Oh")
    document["elements"]["classes"].pop()
    expect(document, "AG-2")


def test_cl1_catches_an_element_moved_to_another_class() -> None:
    document = corrupted("C3v")
    document["elements"]["classes"][1]["members"].append(3)
    document["elements"]["classes"][2]["members"].remove(3)
    expect(document, "CL-1")


def test_cl2_catches_a_bad_or_duplicate_label() -> None:
    document = corrupted("C3v")
    document["elements"]["classes"][1]["label"] = "2C3x"
    expect(document, "CL-2", only=True)
    document["elements"]["classes"][1]["label"] = "E"
    expect(document, "CL-2", only=True)


def test_cl3_catches_a_flipped_regular_flag() -> None:
    for symbol, index in (("Oh", 1), ("C3v", 2), ("D3h", 3)):
        document = corrupted(symbol)
        record = document["elements"]["classes"][index]
        record["regular"] = not record["regular"]
        # DG-6 counts the regular classes
        expect(document, "CL-3", "DG-6", only=True, note=symbol)


def test_cl4_catches_a_wrong_inverse_class() -> None:
    document = corrupted("C6")
    document["elements"]["classes"][1]["inverse"] = 1
    expect(document, "CL-4", only=True)


def test_cl5_catches_swapped_polyhedral_labels() -> None:
    document = corrupted("Td")
    classes = document["elements"]["classes"]
    classes[1]["label"], classes[2]["label"] = classes[2]["label"], classes[1]["label"]
    expect(document, "CL-5")


def test_ag3_catches_a_reordered_parity_block() -> None:
    document = corrupted("C3v")
    document["elements"]["items"][0]["parity"] = -1
    expect(document, "AG-3")


def test_clean_models_pass_every_applicable_check() -> None:
    for symbol in ("C1", "C3v", "D6h", "Oh", "Ih"):
        document = model(symbol)
        assert checks.check(document, applicable(document)) == []


# ---- DG-* and the G-dagger side of CL-* -------------------------------------


def test_dg2_catches_a_flipped_factor_sign() -> None:
    """``MUT-2``: one ``omega(g, h)`` sign flipped is caught by ``DG-2``."""
    for symbol in ("C3v", "Oh"):
        document = corrupted(symbol)
        document["factor"]["table"][1][2] *= -1
        expect(document, "DG-2", note=symbol)
        document = corrupted_double(symbol)
        document["factor"]["table"][1][2] *= -1
        expect(document, "DG-2", note=symbol)


def test_dg1_catches_a_bad_value_and_a_denormalised_table() -> None:
    document = corrupted("C3v")
    document["factor"]["table"][0][1] = -1
    expect(document, "DG-1")
    document = corrupted("C3v")
    document["factor"]["table"][2][3] = 0
    expect(document, "DG-1")


def test_dg3_catches_an_asymmetric_inverse_pair() -> None:
    document = corrupted("C6")
    labels = [record["label"] for record in document["elements"]["items"]]
    c6, c6_inverse = labels.index("C6^1"), labels.index("C6^-1")
    document["factor"]["table"][c6][c6_inverse] *= -1
    expect(document, "DG-3")


def test_dg4_catches_a_corrupted_double_group_table() -> None:
    document = corrupted_double("C3v")
    N = 6
    document["cayley"]["table"][N][1] = 1  # ~E no longer bars C3
    expect(document, "DG-4")
    document = corrupted_double("C3v")
    document["elements"]["items"][7]["label"] = "~C3^-1"
    expect(document, "DG-4")
    document = corrupted_double("C3v")
    document["elements"]["items"][7]["base_index"] = 2
    expect(document, "DG-4")


def test_dg5_catches_a_binary_rotation_that_does_not_square_to_the_full_turn() -> None:
    document = corrupted("D2h")
    labels = [record["label"] for record in document["elements"]["items"]]
    c2 = labels.index("C2")
    document["factor"]["table"][c2][c2] = 1
    expect(document, "DG-5")
    document = corrupted("D2h")
    sigma, c2 = labels.index("sigma_h"), labels.index("C2")
    row = document["factor"]["table"][sigma]
    row[c2] *= -1  # omega(sigma_h, C2) no longer equals omega(C2, C2)
    expect(document, "DG-5")


def test_dg6_catches_a_missing_double_class_and_a_wrong_regular_count() -> None:
    document = corrupted_double("D3h")
    document["elements"]["classes"].pop()
    expect(document, "DG-6")
    document = corrupted("D3h")
    document["elements"]["classes"][2]["regular"] = True
    expect(document, "DG-6")


def test_cl3_catches_a_flipped_split_flag_of_the_double_group() -> None:
    document = corrupted_double("D4")
    classes = document["elements"]["classes"]
    merged = next(cls for cls in classes if not cls["split"])
    merged["split"] = True
    merged["barred"] = False
    expect(document, "CL-3")
    document = corrupted_double("D4")
    split = next(cls for cls in document["elements"]["classes"] if cls["split"])
    split["barred"] = not split["barred"]
    expect(document, "CL-3")


def test_cl4_catches_a_wrong_base_class_of_the_double_group() -> None:
    document = corrupted_double("C3v")
    classes = document["elements"]["classes"]
    classes[2]["base_class"], classes[3]["base_class"] = 0, 0
    expect(document, "CL-4")
    document = corrupted_double("C3v")
    classes = document["elements"]["classes"]
    del classes[3]["barred"]
    expect(document, "CL-4")


def test_clean_double_models_pass_every_applicable_check() -> None:
    for symbol in ("C1", "C3v", "D6h", "Oh", "Ih"):
        document = double_model(symbol)
        assert checks.check(document, applicable(document)) == []


# ---- IL-*, CH-* and VL-* --------------------------------------------------------


def test_il2_catches_swapped_symbols_of_the_polyhedral_rules() -> None:
    """``MUT-5`` on the symbol pairs the polyhedral rules tell apart."""
    document = corrupted("C3v")
    irreps = document["characters"]["irreps"]
    for key in ("symbol", "letter", "subscript"):
        irreps[0][key], irreps[1][key] = irreps[1].get(key), irreps[0].get(key)
    expect(document, "IL-2")
    document = corrupted("D3")
    irreps = document["characters"]["irreps"]
    irreps[4]["symbol"], irreps[5]["symbol"] = irreps[5]["symbol"], irreps[4]["symbol"]
    irreps[4]["component"], irreps[5]["component"] = 2, 1
    irreps[4]["partner"], irreps[5]["partner"] = 5, 4
    expect(document, "IL-2", "IL-4")
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
        expect(document, "IL-2", note=symbol)


def test_il1_catches_a_bad_duplicate_or_inconsistent_symbol() -> None:
    document = corrupted("C3v")
    document["characters"]["irreps"][0]["symbol"] = "F1"
    expect(document, "IL-1")
    document = corrupted("C3v")
    document["characters"]["irreps"][0]["symbol"] = "A2"
    expect(document, "IL-1")
    document = corrupted("C3v")
    document["characters"]["irreps"][0]["letter"] = "B"
    expect(document, "IL-1")


def test_il3_catches_a_wrong_dimension_indicator_or_partner() -> None:
    document = corrupted("C3v")
    document["characters"]["irreps"][2]["dimension"] = 1
    expect(document, "IL-3")
    document = corrupted("C3v")
    document["characters"]["irreps"][3]["frobenius_schur"] = 1
    expect(document, "IL-3", "CH-5")
    document = corrupted("C3v")
    document["characters"]["irreps"][4]["partner"] = 4
    expect(document, "IL-3")
    document = corrupted("C3v")
    del document["characters"]["irreps"][4]["partner"]
    del document["characters"]["irreps"][4]["component"]
    expect(document, "IL-3")


def test_ch3_catches_a_missing_row() -> None:
    document = corrupted("C3v")
    table = document["characters"]
    for key in ("irreps", "exact", "numeric", "pretty"):
        table[key].pop()
    expect(document, "CH-3")


def test_ch4_and_ch6_catch_a_projective_value_on_a_merged_class() -> None:
    document = corrupted("D4")
    table = document["characters"]
    merged = next(c for c, cls in enumerate(table["classes"]) if not cls["split"])
    projective = next(k for k, r in enumerate(table["irreps"]) if r["projective"])
    table["exact"][projective][merged] = table["exact"][0][0]
    expect(document, "CH-4", "CH-6")


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
    expect(document, "CH-4")


def test_vl1_and_vl2_catch_bad_value_forms() -> None:
    document = corrupted("C3v")
    document["characters"]["exact"][0][0] = [[2, 2], [0, 1], [0, 1], [0, 1]]
    expect(document, "VL-1")
    document = corrupted("C3v")
    document["characters"]["exact"][0][0] = [[1, 1], [0, 1]]
    expect(document, "VL-1")
    document = corrupted("C3v")
    document["characters"]["numeric"][0][0] = [1.0000000000000002, 0.0]
    expect(document, "VL-2")
    document = corrupted("C3v")
    document["characters"]["numeric"][0][0] = [1.0, -0.0]
    expect(document, "VL-2")
    document = corrupted("C3v")
    document["characters"]["pretty"][0][0] = "2"
    expect(document, "VL-2")
    document = corrupted("C3v")
    document["characters"]["pretty"][0][0] = "sqrt2"
    expect(document, "VL-2")


def test_checks_report_a_characters_table_without_an_exact_form() -> None:
    document = corrupted("C3v")
    del document["characters"]["exact"]
    expect(document, "IL-2", "CH-1")


# ---- the matrices table: MUT-2, MUT-3 and the other matrix checkers -----------------


def test_mut3_a_matrix_multiplied_by_a_phase_is_caught_by_mx1_or_mx2_and_mx4() -> None:
    document = corrupted("C3v")
    F = cyclotomic_field(document["field"]["conductor"])
    entry = document["matrices"]["irreps"][2]  # E
    entry["exact"][1] = [
        [F.exact(F.mul(F.zeta(1), F.from_exact(pairs))) for pairs in row]
        for row in entry["exact"][1]
    ]
    expect(document, "MX-4", any_of=("MX-1", "MX-2"))


def test_mx3_catches_a_matrix_that_is_not_unitary() -> None:
    document = corrupted("C3v")
    entry = document["matrices"]["irreps"][2]
    entry["exact"][1] = [
        [[[2 * p, q] for p, q in pairs] for pairs in row] for row in entry["exact"][1]
    ]
    expect(document, "MX-3")


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
    expect(document, "MX-5", only=True)


def test_mx7_catches_a_missing_matrix_and_a_foreign_entry() -> None:
    document = corrupted("C3v")
    document["matrices"]["irreps"][0]["exact"].pop()
    expect(document, "MX-7")
    document = corrupted("C3v")
    document["matrices"]["irreps"][0]["exact"][0][0][0] = [
        [1, 1],
        [0, 1],
    ]  # degree 2, not 4
    expect(document, "MX-7")


def test_mx8_catches_a_duplicated_irrep() -> None:
    document = corrupted("C3v")
    entries = document["matrices"]["irreps"]
    entries[1]["exact"] = copy.deepcopy(entries[0]["exact"])
    entries[1]["numeric"] = copy.deepcopy(entries[0]["numeric"])
    expect(document, "MX-4", "MX-8")


# ---- the ledger's mutations, MUT-1 ... MUT-7, on C3v and Oh -----------------------


def _different_entry(table: list[list[int]], row: int, column: int) -> int:
    return (table[row][column] + 1) % len(table)


@pytest.mark.parametrize("symbol", MUTATION_GROUPS)
def test_mut1_a_replaced_cayley_entry_is_caught_by_ag1(symbol: str) -> None:
    document = corrupted(symbol)
    table = document["cayley"]["table"]
    table[1][2] = _different_entry(table, 1, 2)
    expect(document, "AG-1", note=symbol)


@pytest.mark.parametrize("symbol", MUTATION_GROUPS)
def test_mut2_a_flipped_factor_sign_is_caught_by_dg2_and_mx2(symbol: str) -> None:
    document = corrupted(symbol)
    document["factor"]["table"][1][2] *= -1
    expect(document, "DG-2", "MX-2", note=symbol)


def _phase_multiplied(document: dict[str, Any], irrep: int, element: int) -> None:
    """Multiply one matrix of one irrep by the primitive root of the document's field."""
    F = cyclotomic_field(document["field"]["conductor"])
    entry = document["matrices"]["irreps"][irrep]
    entry["exact"][element] = [
        [F.exact(F.mul(F.zeta(1), F.from_exact(pairs))) for pairs in row]
        for row in entry["exact"][element]
    ]


@pytest.mark.parametrize("symbol", MUTATION_GROUPS)
def test_mut3_a_phase_multiplied_matrix_is_caught_by_mx1_or_mx2_and_mx4(
    symbol: str,
) -> None:
    document = corrupted(symbol)
    _phase_multiplied(document, 2, 1)  # a two-dimensional irrep on a generator
    expect(document, "MX-4", any_of=("MX-1", "MX-2"), note=symbol)


@pytest.mark.parametrize("symbol", MUTATION_GROUPS)
def test_mx1_a_corrupted_non_generator_matrix_is_caught_by_closure(
    symbol: str,
) -> None:
    # the last element of the double group is a product of generators, so two
    # words reach it and MX-1 sees the mismatch before the all-pairs MX-2
    document = corrupted(symbol)
    last = len(document["matrices"]["irreps"][2]["exact"]) - 1
    _phase_multiplied(document, 2, last)
    expect(document, "MX-1", "MX-2", note=symbol)


@pytest.mark.parametrize("symbol", MUTATION_GROUPS)
def test_mut4_a_replaced_character_is_caught_by_ch1_and_ch2(symbol: str) -> None:
    document = corrupted(symbol)
    table = document["characters"]
    row, column = next(
        (k, c)
        for k, values in enumerate(table["exact"])
        for c, value in enumerate(values)
        if value != values[0]
    )
    table["exact"][row][column] = table["exact"][row][0]  # the dimension: in the field
    expect(document, "CH-1", "CH-2", note=symbol)


@pytest.mark.parametrize("symbol", MUTATION_GROUPS)
def test_mut5_swapped_class_labels_are_caught_by_cl2_or_cl4(symbol: str) -> None:
    document = corrupted(symbol)
    classes = document["elements"]["classes"]
    classes[1]["label"], classes[2]["label"] = classes[2]["label"], classes[1]["label"]
    expect(document, any_of=("CL-2", "CL-4"), note=symbol)
    document = corrupted_double(symbol)
    classes = document["elements"]["classes"]
    classes[2]["label"], classes[4]["label"] = classes[4]["label"], classes[2]["label"]
    expect(document, any_of=("CL-2", "CL-4"), note=symbol)


@pytest.mark.parametrize("symbol", MUTATION_GROUPS)
def test_mut5_swapped_irrep_symbols_are_caught_by_il2(symbol: str) -> None:
    document = corrupted(symbol)
    first, second = document["characters"]["irreps"][:2]
    for key in set(first) | set(second):
        if key in ("index", "dimension", "frobenius_schur", "projective"):
            continue
        a, b = first.pop(key, None), second.pop(key, None)
        if b is not None:
            first[key] = b
        if a is not None:
            second[key] = a
    expect(document, "IL-2", note=symbol)


@pytest.mark.parametrize("symbol", MUTATION_GROUPS)
def test_mut6_a_moved_element_is_caught_by_cl1_and_ch4(symbol: str) -> None:
    document = corrupted(symbol)
    classes = document["elements"]["classes"]
    moved = classes[2]["members"][0]
    classes[2]["members"].remove(moved)
    classes[1]["members"].append(moved)
    document["elements"]["items"][moved]["class"] = 1
    expect(document, "CL-1", "CH-4", note=symbol)
    document = corrupted_double(symbol)
    classes = document["elements"]["classes"]
    moved = classes[2]["members"][0]
    classes[2]["members"].remove(moved)
    classes[1]["members"].append(moved)
    document["elements"]["items"][moved]["class"] = 1
    expect(document, "CL-1", "CH-4", note=symbol)


@pytest.mark.reference
@pytest.mark.parametrize("symbol", MUTATION_GROUPS)
def test_mut7_an_altered_reference_row_is_caught_by_ref1(symbol: str) -> None:
    """The structural comparison of ``REF-1`` finds no bijection once one row of
    the transcribed Altmann & Herzig table is altered."""
    import math
    import tomllib
    from pathlib import Path

    from test_reference_characters import Computed, structural_match

    from pointgroup_tables.values import parse

    path = (
        Path(__file__).parent / "reference/altmann-herzig/characters" / f"{symbol}.toml"
    )
    if not path.is_file():
        pytest.skip(f"missing reference file {path}")
    data = tomllib.loads(path.read_text())
    comp = Computed(symbol, double=False)
    F_book = cyclotomic_field(data["conductor"])
    F = cyclotomic_field(math.lcm(F_book.conductor, comp.F.conductor))
    book_classes = [(c["size"], c["order"], c["parity"]) for c in data["classes"]]
    labels = [c["label"] for c in data["classes"]]
    book_rows = [
        (
            row["label"],
            [F_book.lift(parse(F_book, text), F) for text in row["characters"]],
            row["projective"],
        )
        for row in data["irreps"]
    ]
    for entry in data.get("known_discrepancy", []):
        a, b = (labels.index(label) for label in entry["swap_projective"])
        for label, values, projective in book_rows:
            if projective:
                values[a], values[b] = values[b], values[a]
    comp_rows = [
        (text, [comp.F.lift(x, F) for x in values], projective)
        for text, values, projective, _ in comp.rows
    ]

    def matched() -> bool:
        found = structural_match(
            book_classes, book_rows, comp.invariants, comp_rows, F, comp.regular
        )
        return found is not None

    assert matched()
    _, values, _ = book_rows[-1]
    values[1] = F.neg(values[1]) if not F.is_zero(values[1]) else values[0]
    assert not matched()


# ---- every check --verify runs is the target of a mutation ------------------------


def test_every_verify_check_is_the_target_of_a_mutation(
    request: pytest.FixtureRequest,
) -> None:
    """Runs last: with the whole module collected, the targets recorded by
    ``expect`` cover every check of the ``--verify`` column of the ledger."""
    import sys
    from pathlib import Path

    collected = {
        getattr(item, "originalname", item.name)
        for item in request.session.items
        if getattr(item, "path", None) == Path(__file__)
    }
    defined = {name for name in dir(sys.modules[__name__]) if name.startswith("test_")}
    if collected != defined:
        pytest.skip("the coverage of the mutations needs the whole module collected")
    missing = set(checks.VERIFY_CHECKS) - TARGETED
    assert not missing, f"no mutation targets {sorted(missing)}"
