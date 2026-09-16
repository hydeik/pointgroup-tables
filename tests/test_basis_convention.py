"""The ``altmann-herzig`` basis convention applied to a document
(docs/spec/basis-conventions.md "Reference basis conventions", "The
spin-harmonic rule", "Verification"; docs/spec/output-schema.md, coverage):
the listing, the metadata and the field, the ``MX-*`` checks under the
convention on ``default``, the intertwiner with the canonical basis on
``seed``, the closed-form examples of the chapter, the polyhedral generator
images, and the uncovered-group error with a fixture convention.

These checks use the shipped convention file only; the comparisons with the
transcriptions are ``REF-7`` in ``tests/test_reference_basis.py``.
"""

from __future__ import annotations

import functools
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

import jsonschema
import pytest
from conftest import DEFAULT_GROUPS, SEED_GROUPS, default_and_large
from test_elements import SCHEMA
from test_matrices import _assert_union, _mx2_groups

from pointgroup_tables import basis, checks, conventions
from pointgroup_tables.abstract import point_group
from pointgroup_tables.cli import main
from pointgroup_tables.document import (
    Request,
    build_document,
    build_model,
    model_document,
)
from pointgroup_tables.field import Matrix, cyclotomic_field
from pointgroup_tables.irreps import named_elements
from pointgroup_tables.listing import conventions_document
from pointgroup_tables.symbols import parse_group
from pointgroup_tables.verify import verification_checks, verify

AH = "altmann-herzig"
# the matrix checks that run plainly under the convention (MX-5 is canonical-only;
# VL-3, the minimal conductor, is a property of the canonical basis: the book's
# 2T matrices are declared in Q(zeta_24) though they lie in Q(zeta_12))
CHECKS = ["MX-1", "MX-3", "MX-4", "MX-6", "MX-7", "MX-8"]


@functools.cache
def model(symbol: str) -> dict[str, Any]:
    """The full document model of a group under the convention, built once."""
    return model_document(symbol, basis_convention=AH)


# ---- the listing and the metadata --------------------------------------------------------


def test_the_shipped_basis_conventions_are_listed() -> None:
    assert list(conventions.BASIS_CONVENTIONS) == ["canonical", AH]
    block = conventions.basis_convention(AH)
    assert list(block) == ["name", "source"]
    assert "Altmann" in block["source"]
    listed = conventions_document()["basis_conventions"]
    assert listed[0] == {"name": "canonical"}
    assert listed[1] == block


def test_the_rules_of_the_shipped_conventions() -> None:
    canonical = conventions.basis_convention_rules("canonical")
    assert not canonical.reference
    assert canonical.axial_rule == "canonical" and canonical.anchors is None
    assert canonical.covers == {"2T", "2O", "2I"}
    ah = conventions.basis_convention_rules(AH)
    assert ah.reference
    assert ah.axial_rule == "spin-harmonic" and ah.anchors == AH
    assert ah.covers == {"2T", "2O", "2I"}
    assert (ah.axial_max_n, ah.polyhedral) == (10, True)
    for cover, conductor, generators in (
        ("2T", 24, ("C3^1(1)", "C2(1)")),
        ("2O", 24, ("C4^1(1)", "~C2'(1)")),
        ("2I", 120, ("C5^1(1)", "C3^1(7)", "C2(4)")),
    ):
        data = conventions.polyhedral_data(AH, cover)
        assert data.conductor == conductor
        assert data.generators == generators
        assert all(len(images) == len(generators) for images in data.irreps.values())


@pytest.mark.parametrize(
    ("symbol", "coverage"),
    [
        ("D6h", "reference"),
        ("Ih", "reference"),
        ("C3v", "reference"),
        ("C10v", "reference"),
        ("D14h", "extrapolated"),
        ("S22", "extrapolated"),
    ],
)
def test_coverage_follows_the_convention_range(symbol: str, coverage: str) -> None:
    document = build_document(Request("matrices", symbol, basis_convention=AH))
    assert list(document["basis_convention"]) == ["name", "source", "coverage"]
    assert document["basis_convention"]["coverage"] == coverage
    assert list(document)[:8] == [
        "schema_version", "generator", "group", "tables", "double", "irreps",
        "values", "field",
    ]  # fmt: skip
    assert list(document)[8:] == ["basis_convention", "elements", "matrices"]


def test_both_conventions_are_reported_side_by_side() -> None:
    document = build_document(
        Request("matrices", "D6h", label_convention=AH, basis_convention=AH)
    )
    assert document["label_convention"]["coverage"] == "reference"
    assert document["basis_convention"]["coverage"] == "reference"
    assert list(document)[8:10] == ["label_convention", "basis_convention"]
    assert model_document("D6h")["basis_convention"] == {"name": "canonical"}


@pytest.mark.parametrize(
    ("symbol", "conductor"),
    [
        ("T", 24), ("Td", 24), ("Th", 24), ("O", 24), ("Oh", 24), ("I", 120), ("Ih", 120),
        ("C3v", 12), ("D2d", 8), ("D3h", 12), ("D6h", 12), ("C5v", 20), ("D8", 16),
        ("C7", 7), ("S8", 16), ("D9", 36), ("D14h", 28),
    ],
)  # fmt: skip
def test_the_field_is_the_declared_or_the_canonical_conductor(
    symbol: str, conductor: int
) -> None:
    document = build_document(Request("matrices", symbol, basis_convention=AH))
    assert document["field"]["conductor"] == conductor
    if parse_group(symbol).polyhedral:
        assert conductor != model_document(symbol)["field"]["conductor"]
    else:
        assert conductor == model_document(symbol)["field"]["conductor"]


# ---- the checks under the convention -------------------------------------------------------


@pytest.mark.parametrize("check_id", CHECKS)
@pytest.mark.parametrize("symbol", default_and_large())
def test_matrix_checks_under_the_convention(check_id: str, symbol: str) -> None:
    assert checks.check(model(symbol), [check_id]) == []


@pytest.mark.parametrize("symbol", _mx2_groups())
def test_mx2_under_the_convention(symbol: str) -> None:
    assert checks.check(model(symbol), ["MX-2"]) == []


@pytest.mark.parametrize("symbol", ["C3v", "D2d", "D6h", "T", "O"])
def test_every_applicable_check_passes_under_the_convention(symbol: str) -> None:
    document = model(symbol)
    names = [name for name in checks.applicable_checks(document) if name != "VL-3"]
    assert checks.check(document, names) == []
    assert checks.check(document, ["MX-5"]) == []  # canonical-only: vacuous here


def test_mx6_is_vacuous_on_the_canonical_basis() -> None:
    assert checks.check(model_document("O"), ["MX-6"]) == []
    assert "MX-6" in verification_checks(Request("matrices", "O"))
    assert "MX-6" in checks.VERIFY_CHECKS


# ---- the intertwiner with the canonical basis ----------------------------------------------


def _is_monomial(F: Any, U: Matrix) -> bool:
    """One non-zero entry per row and per column."""
    d = len(U)
    rows = [sum(not F.is_zero(x) for x in row) for row in U]
    columns = [sum(not F.is_zero(U[i][j]) for i in range(d)) for j in range(d)]
    return rows == [1] * d and columns == [1] * d


@pytest.mark.parametrize("symbol", SEED_GROUPS)
def test_the_intertwiner_exists_and_is_unitary(symbol: str) -> None:
    """Computed here from the two realisations, never stored: ``U D_can(g) = D(g) U``
    with ``U^dagger U`` a non-zero real scalar (positive, being ``U^dagger U``); a
    permutation times a diagonal phase for the two-dimensional axial irreps (the
    consequence the chapter states)."""
    spec = parse_group(symbol)
    reference, canonical = build_model(spec, AH), build_model(spec)
    F, F_can = reference.field, canonical.field
    by_symbol = {irrep.symbol.text: irrep for irrep in canonical.irreps}
    for irrep in reference.irreps:
        if irrep.dimension < 2:
            continue
        lifted = [
            [[F_can.lift(x, F) for x in row] for row in matrix]
            for matrix in by_symbol[irrep.symbol.text].matrices
        ]
        U = basis.intertwiner(F, irrep.matrices, lifted)
        assert U is not None, irrep.symbol.text
        scale = basis.unitary_scale(F, U)
        assert scale is not None and F.is_real(scale), irrep.symbol.text
        if not spec.polyhedral:
            assert _is_monomial(F, U), irrep.symbol.text


def test_the_intertwiner_of_the_canonical_basis_with_itself_is_scalar() -> None:
    canonical = build_model(parse_group("D3h"))
    F = canonical.field
    for irrep in canonical.irreps:
        if irrep.dimension < 2:
            continue
        U = basis.intertwiner(F, irrep.matrices, irrep.matrices)
        assert U is not None and basis.unitary_scale(F, U) is not None
        assert (
            all(F.is_zero(U[i][j]) for i in range(2) for j in range(2) if i != j)
            and U[0][0] == U[1][1]
        )


def test_inequivalent_realisations_have_no_intertwiner() -> None:
    canonical = build_model(parse_group("D3h"))
    F = canonical.field
    two = [irrep for irrep in canonical.irreps if irrep.dimension == 2]
    assert basis.intertwiner(F, two[0].matrices, two[1].matrices) is None


# ---- the spin-harmonic rule --------------------------------------------------------------


@pytest.mark.parametrize(
    ("N", "levels"),
    [(2, [2]), (3, [3]), (4, [2, 4]), (6, [3, 6]), (8, [2, 4, 8]), (9, [3, 9]),
     (10, [5, 10]), (12, [3, 6, 12]), (15, [5, 15]), (18, [3, 9, 18]), (24, [3, 6, 12, 24])],
)  # fmt: skip
def test_the_chain_divides_by_the_smallest_prime_factor(
    N: int, levels: list[int]
) -> None:
    assert basis.chain(N) == levels


@pytest.mark.parametrize(
    ("N", "k2", "m1", "residue"),
    [
        # the examples of the chapter (D6: E1, E2, E5/2, E3/2; D9: E5/2, E3)
        (6, 2, Fraction(1), Fraction(1)),
        (6, 4, Fraction(-2), Fraction(1)),
        (6, 5, Fraction(-5, 2), Fraction(1, 2)),
        (6, 3, Fraction(3, 2), Fraction(3, 2)),
        (9, 5, Fraction(-5, 2), Fraction(1, 2)),
        (9, 6, Fraction(3), Fraction(3)),
        # E_{N/4}: -k for N/4 odd and greater than one, +k otherwise
        (12, 6, Fraction(-3), Fraction(-3)),
        (20, 10, Fraction(-5), Fraction(-5)),
        (8, 4, Fraction(2), Fraction(2)),
        (16, 8, Fraction(4), Fraction(4)),
        (4, 2, Fraction(1), Fraction(1)),
    ],
)
def test_the_first_partner_and_the_deciding_residue(
    N: int, k2: int, m1: Fraction, residue: Fraction
) -> None:
    m2, mt2 = basis.first_partner(N, k2)
    assert (Fraction(m2, 2), Fraction(mt2, 2)) == (m1, residue)


def _images(symbol: str) -> tuple[Any, dict[str, list[Matrix]], int, int]:
    """``(F, symbol -> matrices, a, b)``: the matrices of every irrep of the group under
    the convention and the indices of the principal and the secondary element."""
    spec = parse_group(symbol)
    reference = build_model(spec, AH)
    named = named_elements(point_group(spec))
    assert named.secondary is not None
    return (
        reference.field,
        {irrep.symbol.text: irrep.matrices for irrep in reference.irreps},
        named.principal,
        named.secondary,
    )


def test_the_examples_of_the_chapter_in_d6() -> None:
    F, D, a, b = _images("D6")
    u = F.zeta  # zeta_12
    zero, one, i = F.zero, F.one, F.constant("i")
    swap = [[zero, one], [one, zero]]
    assert D["E1"][a] == [[u(-2 % 12), zero], [zero, u(2)]]
    assert D["E1"][b] == [[zero, F.neg(one)], [F.neg(one), zero]]
    assert D["E2"][a] == [[u(4), zero], [zero, u(-4 % 12)]]
    assert D["E2"][b] == D["E1"][b]
    assert D["E5/2"][a] == [[u(5), zero], [zero, u(-5 % 12)]]
    assert D["E5/2"][b] == [[F.neg(x) for x in row] for row in [[zero, i], [i, zero]]]
    assert D["E3/2"][b] == [[zero, i], [i, zero]]
    assert D["A1"][b] == [[one]] and D["A2"][b] == [[F.neg(one)]]
    assert swap != D["E1"][b]


def test_the_examples_of_the_chapter_in_c6v_d9_and_c5v() -> None:
    F, D, _, b = _images("C6v")
    zero, one = F.zero, F.one
    assert D["E1"][b] == [[zero, one], [one, zero]]
    assert D["E1/2"][b] == [[zero, F.neg(one)], [one, zero]]
    F, D, a, b = _images("D9")
    u = F.zeta  # zeta_36 here; u = zeta_18 is zeta_36^2
    zero, one = F.zero, F.one
    assert D["E5/2"][a] == [[u(10), zero], [zero, u(-10 % 36)]]
    assert D["E3"][a] == [[u(-12 % 36), zero], [zero, u(12)]]
    assert D["E3"][b] == [[zero, F.neg(one)], [F.neg(one), zero]]
    F, D, _, b = _images("C5v")
    zero, one = F.zero, F.one
    assert D["E1"][b] == [[zero, F.neg(one)], [F.neg(one), zero]]
    assert D["E1/2"][b] == [[zero, F.neg(one)], [one, zero]]


def test_the_secondary_azimuth_follows_the_anchoring() -> None:
    quarter = {
        symbol: basis.secondary_azimuth(AH, parse_group(symbol))
        for symbol in ("C2v", "C3v", "C4v", "C5v", "C6v", "C8v", "C10v", "C12v",
                       "C14v", "D6", "D6h", "D6d", "D3h", "D2d", "D14h")
    }  # fmt: skip
    assert {symbol for symbol, q in quarter.items() if q == 1} == {
        "C2v",
        "C3v",
        "C5v",
        "C6v",
        "C10v",
        "C14v",
    }
    with pytest.raises(Exception, match="unknown label convention"):
        basis.secondary_azimuth("bilbao", parse_group("C6v"))


def test_cyclic_rotation_groups_are_unchanged_by_the_rule() -> None:
    for symbol in ("C6", "C6h", "S8", "C7"):
        ours = build_document(Request("matrices", symbol, basis_convention=AH))
        canonical = build_document(Request("matrices", symbol))
        assert ours["matrices"] == canonical["matrices"]


def test_the_canonical_rule_through_the_same_functions() -> None:
    """``axial_images`` of the canonical convention are the closed forms of the
    canonical basis (the images ``test_matrices`` checks on ``C3v``)."""
    G = point_group(parse_group("D3"))
    F = cyclotomic_field(12)
    Da, Db = basis.axial_images("canonical", G, F, 2)
    assert Da == [[F.zeta(4), F.zero], [F.zero, F.zeta(8)]]
    assert Db == [[F.zero, F.one], [F.one, F.zero]]
    Da, Db = basis.axial_images("canonical", G, F, 1)
    i = F.constant("i")
    assert Da == [[F.zeta(2), F.zero], [F.zero, F.zeta(10)]]
    assert Db == [[F.zero, i], [i, F.zero]]


# ---- the polyhedral data ---------------------------------------------------------------


def test_the_book_s_generator_images_are_kept_verbatim() -> None:
    """The lifted generators receive the file's matrices as they stand (the
    barred lift of ``C2'(1)`` of ``2O`` included: there the document's matrix of
    the element is ``-D`` in the projective irreps), the other elements follow by
    closure."""
    for symbol, cover, R in (("T", "2T", "T"), ("O", "2O", "O"), ("I", "2I", "I")):
        reference = build_model(parse_group(symbol), AH)
        F = reference.field
        data = conventions.polyhedral_data(AH, cover)
        generators, _ = basis.lifted_generators(AH, reference.double.cover, R)
        lifted = [q for q, _ in reference.double.elements[: reference.order]]
        ebar = reference.double.cover.ebar
        for irrep in reference.irreps:
            if irrep.symbol.text not in data.irreps:
                continue
            for q, image in zip(
                generators, data.irreps[irrep.symbol.text], strict=True
            ):
                barred = q not in lifted
                g = lifted.index(reference.double.cover.mul(ebar, q) if barred else q)
                expected = image
                if barred and irrep.projective:
                    expected = [[F.neg(x) for x in row] for row in image]
                assert irrep.matrices[g] == expected, (symbol, irrep.symbol.text)
    O = conventions.polyhedral_data(AH, "2O")
    assert O.generators[1] == "~C2'(1)"


def test_type_ii_and_twisted_groups_share_the_data_of_the_rotation_group() -> None:
    """``Td`` (twisted) and ``Oh`` (centrosymmetric) realise every irrep by the
    matrices of ``O`` on the rotation part, with the parity sign for ``Oh``."""
    O = build_model(parse_group("O"), AH)
    index_of = {r: k for k, (r, _) in enumerate(O.group.elements)}
    by_symbol_O = {irrep.symbol.text: irrep for irrep in O.irreps}
    for symbol in ("Td", "Oh"):
        model_ = build_model(parse_group(symbol), AH)
        assert model_.field == O.field
        for irrep in model_.irreps:
            key = irrep.symbol.text
            sign = 1
            if symbol == "Oh":
                key, sign = key[:-1], (1 if key.endswith("g") else -1)
            source = by_symbol_O[key]
            for k, (r, parity) in enumerate(model_.group.elements):
                expected = source.matrices[index_of[r]]
                if symbol == "Oh" and parity == -1 and sign == -1:
                    expected = [[O.field.neg(x) for x in row] for row in expected]
                assert irrep.matrices[k] == expected, (symbol, key)


# ---- the uncovered-group error ----------------------------------------------------------


@pytest.fixture
def fixture_convention(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> str:
    """A basis convention with the canonical polyhedral data and no ``[axial]`` table."""
    name = "fixture-polyhedral-only"
    text = conventions.BASIS_FILES["canonical"].read_text(encoding="utf-8")
    text = text.replace('name = "canonical"', f'name = "{name}"\nsource = "a fixture"')
    text = text.replace('[axial]\nrule = "canonical"\n', "")
    assert "[axial]" not in text
    path = tmp_path / f"{name}.toml"
    path.write_text(text, encoding="utf-8")
    monkeypatch.setitem(conventions.BASIS_FILES, name, path)
    monkeypatch.setitem(
        conventions.BASIS_CONVENTIONS, name, {"name": name, "source": "a fixture"}
    )
    return name


def test_a_convention_without_axial_coverage_exits_4(
    fixture_convention: str, tmp_path: Path, capsysbinary: pytest.CaptureFixture[bytes]
) -> None:
    rules = conventions.basis_convention_rules(fixture_convention)
    assert rules.axial_rule is None and rules.covers == {"2T", "2O", "2I"}
    out = tmp_path / "C3v.json"
    assert (
        main(
            [
                "matrices",
                "C3v",
                "--basis-convention",
                fixture_convention,
                "-o",
                str(out),
            ]
        )
        == 4
    )
    captured = capsysbinary.readouterr()
    assert captured.out == b"" and not out.exists()
    assert captured.err.decode() == (
        f"pointgroup-tables: error: basis convention '{fixture_convention}' does not cover C3v\n"
        "hint: run 'pointgroup-tables list conventions' for the shipped conventions\n"
    )
    assert main(["elements", "D6h", "--basis-convention", fixture_convention]) == 4
    assert main(["matrices", "Oh", "--basis-convention", fixture_convention]) == 0
    document = json.loads(capsysbinary.readouterr().out)
    assert document["basis_convention"] == {
        "name": fixture_convention,
        "source": "a fixture",
        "coverage": "extrapolated",
    }
    assert document["field"]["conductor"] == 8


def test_an_unknown_basis_convention_exits_4(
    capsysbinary: pytest.CaptureFixture[bytes],
) -> None:
    assert main(["matrices", "C3v", "--basis-convention", "bilbao"]) == 4
    captured = capsysbinary.readouterr()
    assert captured.out == b""
    assert b"unknown basis convention 'bilbao'" in captured.err


# ---- the command line, the schema and the union of all ------------------------------------


def test_the_convention_on_the_command_line(
    capsysbinary: pytest.CaptureFixture[bytes],
) -> None:
    assert main(["matrices", "D6h", "--basis-convention", AH, "--verify"]) == 0
    document = json.loads(capsysbinary.readouterr().out)
    assert document["basis_convention"]["coverage"] == "reference"
    jsonschema.Draft202012Validator(SCHEMA).validate(document)
    assert main(["list", "conventions"]) == 0
    listed = json.loads(capsysbinary.readouterr().out)["basis_conventions"]
    assert [entry["name"] for entry in listed] == ["canonical", AH]


@pytest.mark.parametrize("symbol", ["C6v", "D3h", "D2d", "T", "O"])
def test_verify_passes_under_the_convention(symbol: str) -> None:
    assert verify(Request("matrices", symbol, basis_convention=AH)) == []


@pytest.mark.parametrize("symbol", ["C6v", "D6h", "O", "Th"])
def test_os1_documents_under_the_convention_validate(symbol: str) -> None:
    validator = jsonschema.Draft202012Validator(SCHEMA)
    for subcommand in ("matrices", "characters", "all"):
        validator.validate(
            build_document(Request(subcommand, symbol, basis_convention=AH))
        )


def test_all_is_the_union_under_the_convention() -> None:
    _assert_union(Request("all", "D3h", basis_convention=AH), "json")
    _assert_union(
        Request("all", "Td", basis_convention=AH, label_convention=AH), "toml"
    )


def test_characters_do_not_depend_on_the_basis() -> None:
    """The same irrep records, classes and numeric characters; the same exact
    vectors where the field is the same (a polyhedral group's field grows)."""
    for symbol in DEFAULT_GROUPS[:12] + ["T", "O", "I"]:
        ours = build_document(Request("characters", symbol, basis_convention=AH))
        canonical = build_document(Request("characters", symbol))
        for key in ("class_group", "irreps", "classes", "numeric"):
            assert ours["characters"][key] == canonical["characters"][key], symbol
        if not parse_group(symbol).polyhedral:
            assert ours["characters"] == canonical["characters"]


def test_determinism_under_the_convention(
    capsysbinary: pytest.CaptureFixture[bytes],
) -> None:
    argv = ["all", "D4d", "--basis-convention", AH, "-f", "yaml"]
    assert main(argv) == 0
    first = capsysbinary.readouterr().out
    assert main(argv) == 0
    assert capsysbinary.readouterr().out == first
