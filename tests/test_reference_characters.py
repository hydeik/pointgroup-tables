"""``REF-1`` and ``REF-2``: structural comparison of the computed character tables
with the Altmann & Herzig and Bilbao transcriptions (docs/spec/verification.md,
"Structural comparison").

A structural comparison uses no label: the reference classes carry their
invariants (size, element order, parity), and a bijection of classes
preserving them together with a bijection of irreps must reproduce the table
entry by entry, exactly. The Altmann & Herzig tables list the classes of ``G``
and give the projective rows on the book's untilded elements, so those rows
match up to one sign per regular class; the Bilbao double-group tables list
the classes of ``G-dagger`` and match without any sign freedom.
"""

from __future__ import annotations

import math
import tomllib
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import pytest
from conftest import double_model, model

from pointgroup_tables.document import model_document
from pointgroup_tables.field import CyclotomicField, Value, cyclotomic_field
from pointgroup_tables.values import parse

REFERENCE = Path(__file__).resolve().parent / "reference"
ALTMANN_HERZIG = REFERENCE / "altmann-herzig" / "characters"
BILBAO = REFERENCE / "bilbao" / "characters"
TAU_INDICATOR = {"a": 1, "b": 0, "c": -1}
REALITY_INDICATOR = {1: 1, 0: 0, -1: -1}

Invariants = tuple[int, int, int]


def _files(directory: Path, suffix: str = "") -> list[str]:
    """The groups transcribed under a directory (``-double`` files by ``suffix``)."""
    if not directory.is_dir():
        return []
    return sorted(
        path.stem.removesuffix(suffix)
        for path in directory.glob(f"*{suffix}.toml")
        if suffix or not path.stem.endswith("-double")
    )


# ---- the computed side ----------------------------------------------------------------


def _orders(table: list[list[int]]) -> list[int]:
    orders = []
    for g in range(len(table)):
        n, h = 1, g
        while h != 0:
            h = table[h][g]
            n += 1
        orders.append(n)
    return orders


class Computed:
    """The computed table in the column space of ``G`` (lifted classes) or ``G-dagger``."""

    def __init__(self, symbol: str, double: bool) -> None:
        document = model(symbol)
        self.F = cyclotomic_field(document["field"]["conductor"])
        table = document["characters"]
        items = document["elements"]["items"]
        columns = table["classes"]
        N = document["group"]["order"]
        if double:
            orders = _orders(double_model(symbol)["cayley"]["table"])
            chosen = list(range(len(columns)))
            self.invariants: list[Invariants] = [
                (
                    cls["size"],
                    orders[cls["representative"]],
                    items[cls["representative"] % N]["parity"],
                )
                for cls in columns
            ]
        else:
            orders = _orders(document["cayley"]["table"])
            chosen = [
                c for c, cls in enumerate(columns) if not cls.get("barred", False)
            ]
            self.invariants = [
                (
                    columns[c]["size"]
                    if columns[c]["split"]
                    else columns[c]["size"] // 2,
                    orders[columns[c]["representative"]],
                    items[columns[c]["representative"]]["parity"],
                )
                for c in chosen
            ]
        self.labels = [columns[c]["label"] for c in chosen]
        self.regular = [columns[c]["split"] for c in chosen]
        self.rows: list[tuple[str, list[Value], bool, int]] = [
            (
                record["symbol"],
                [self.F.from_exact(table["exact"][k][c]) for c in chosen],
                record["projective"],
                record["frobenius_schur"],
            )
            for k, record in enumerate(table["irreps"])
        ]


# ---- the matching ------------------------------------------------------------------------


def structural_match(
    book_classes: Sequence[Invariants],
    book_rows: Sequence[tuple[str, list[Value], bool]],
    comp_classes: Sequence[Invariants],
    comp_rows: Sequence[tuple[str, list[Value], bool]],
    F: CyclotomicField,
    sign_freedom: Sequence[bool] | None,
) -> tuple[list[int], list[int], dict[int, int]] | None:
    """``(perm, pairing, signs)``: ``perm[i]`` the computed column of book column ``i``,
    ``pairing[r]`` the computed row of book row ``r``, ``signs`` the sign applied to
    projective rows per computed column (Altmann & Herzig only), or ``None``."""
    n = len(book_classes)
    if sorted(book_classes) != sorted(comp_classes):
        return None
    candidates = [
        [j for j in range(n) if comp_classes[j] == book_classes[i]] for i in range(n)
    ]
    order = sorted(range(n), key=lambda i: len(candidates[i]))

    def row_signs(
        book: list[Value],
        comp: list[Value],
        projective: bool,
        perm: list[int | None],
        signs: dict[int, int],
    ) -> dict[int, int] | None:
        """The signs a pairing of two rows needs on the assigned columns, or ``None``."""
        out = dict(signs)
        for i, j in enumerate(perm):
            if j is None:
                continue
            value, target = book[i], comp[j]
            free = projective and sign_freedom is not None and sign_freedom[j]
            if not free:
                if value != target:
                    return None
                continue
            if F.is_zero(value):
                if not F.is_zero(target):
                    return None
                continue
            if value == target:
                needed = 1
            elif value == F.neg(target):
                needed = -1
            else:
                return None
            if out.get(j, needed) != needed:
                return None
            out[j] = needed
        return out

    def pair_rows(
        r: int, perm: list[int | None], used: set[int], signs: dict[int, int]
    ) -> tuple[list[int], dict[int, int]] | None:
        if r == len(book_rows):
            return [], signs
        _, values, projective = book_rows[r]
        for s, (_, comp_values, comp_projective) in enumerate(comp_rows):
            if s in used or comp_projective != projective:
                continue
            new_signs = row_signs(values, comp_values, projective, perm, signs)
            if new_signs is None:
                continue
            rest = pair_rows(r + 1, perm, used | {s}, new_signs)
            if rest is not None:
                return [s, *rest[0]], rest[1]
        return None

    def search(
        position: int, perm: list[int | None], used: set[int]
    ) -> tuple[list[int], list[int], dict[int, int]] | None:
        if position == n:
            paired = pair_rows(0, perm, set(), {})
            if paired is None:
                return None
            return [j for j in perm if j is not None], paired[0], paired[1]
        i = order[position]
        for j in candidates[i]:
            if j in used:
                continue
            perm[i] = j
            if pair_rows(0, perm, set(), {}) is not None:
                found = search(position + 1, perm, used | {j})
                if found is not None:
                    return found
            perm[i] = None
        return None

    return search(0, [None] * n, set())


def _lift_all(
    values: list[Value], source: CyclotomicField, target: CyclotomicField
) -> list[Value]:
    return [source.lift(x, target) for x in values]


# ---- REF-1: Altmann & Herzig -----------------------------------------------------------------


def _read_altmann_herzig(symbol: str) -> dict[str, Any]:
    path = ALTMANN_HERZIG / f"{symbol}.toml"
    if not path.is_file():
        pytest.skip(f"missing reference file {path}")
    return tomllib.loads(path.read_text())


def _book_rows(
    data: dict[str, Any], F_book: CyclotomicField, F: CyclotomicField
) -> list[tuple[str, list[Value], bool]]:
    return [
        (
            row["label"],
            _lift_all([parse(F_book, text) for text in row["characters"]], F_book, F),
            row["projective"],
        )
        for row in data["irreps"]
    ]


def _swap_projective(
    rows: list[tuple[str, list[Value], bool]], classes: list[str], a: str, b: str
) -> list[tuple[str, list[Value], bool]]:
    ia, ib = classes.index(a), classes.index(b)
    out = []
    for label, values, projective in rows:
        values = list(values)
        if projective:
            values[ia], values[ib] = values[ib], values[ia]
        out.append((label, values, projective))
    return out


@pytest.mark.reference
@pytest.mark.parametrize("symbol", _files(ALTMANN_HERZIG))
def test_ref1_altmann_herzig_character_tables(symbol: str) -> None:
    data = _read_altmann_herzig(symbol)
    comp = Computed(symbol, double=False)
    F_book = cyclotomic_field(data["conductor"])
    F = cyclotomic_field(math.lcm(F_book.conductor, comp.F.conductor))
    book_classes = [(c["size"], c["order"], c["parity"]) for c in data["classes"]]
    book_rows = _book_rows(data, F_book, F)
    comp_rows = [
        (symbol_text, _lift_all(values, comp.F, F), projective)
        for symbol_text, values, projective, _ in comp.rows
    ]
    found = structural_match(
        book_classes, book_rows, comp.invariants, comp_rows, F, comp.regular
    )
    discrepancies = data.get("known_discrepancy", [])
    if discrepancies:
        assert found is None, (
            "the documented discrepancy is not a discrepancy: the printed table matches"
        )
        labels = [c["label"] for c in data["classes"]]
        for entry in discrepancies:
            book_rows = _swap_projective(book_rows, labels, *entry["swap_projective"])
        found = structural_match(
            book_classes, book_rows, comp.invariants, comp_rows, F, comp.regular
        )
    assert found is not None, (
        f"no class and irrep bijection reproduces the book's table of {symbol}"
    )
    _, pairing, _ = found
    assert sorted(pairing) == list(range(len(comp_rows))), (
        "the book lists a different number of irreps"
    )
    for row, s in zip(data["irreps"], pairing, strict=True):
        assert TAU_INDICATOR[row["tau"]] == comp.rows[s][3], (
            row["label"],
            comp.rows[s][0],
        )


# ---- REF-2: Bilbao -----------------------------------------------------------------------------


def _read_bilbao(symbol: str, double: bool) -> dict[str, Any]:
    path = BILBAO / f"{symbol}{'-double' if double else ''}.toml"
    if not path.is_file():
        pytest.skip(f"missing reference file {path}")
    data = tomllib.loads(path.read_text())
    assert data["double"] is double
    return data


@pytest.mark.reference
@pytest.mark.parametrize("symbol", _files(BILBAO))
def test_ref2_bilbao_point_group_tables(symbol: str) -> None:
    data = _read_bilbao(symbol, double=False)
    comp = Computed(symbol, double=False)
    F_book = cyclotomic_field(data["conductor"])
    F = cyclotomic_field(math.lcm(F_book.conductor, comp.F.conductor))
    book_classes = [(c["size"], c["order"], c["parity"]) for c in data["classes"]]
    book_rows = _book_rows(data, F_book, F)
    assert all(not projective for _, _, projective in book_rows)
    comp_rows = [
        (s, _lift_all(values, comp.F, F), projective)
        for s, values, projective, _ in comp.rows
        if not projective
    ]
    found = structural_match(
        book_classes, book_rows, comp.invariants, comp_rows, F, None
    )
    assert found is not None, (
        f"no class and irrep bijection reproduces the server's table of {symbol}"
    )
    assert sorted(found[1]) == list(range(len(comp_rows)))


@pytest.mark.reference
@pytest.mark.parametrize("symbol", _files(BILBAO, "-double"))
def test_ref2_bilbao_double_group_tables(symbol: str) -> None:
    data = _read_bilbao(symbol, double=True)
    comp = Computed(symbol, double=True)
    F_book = cyclotomic_field(data["conductor"])
    F = cyclotomic_field(math.lcm(F_book.conductor, comp.F.conductor))
    book_classes = [(c["size"], c["order"], c["parity"]) for c in data["classes"]]
    book_rows = _book_rows(data, F_book, F)
    comp_rows = [
        (s, _lift_all(values, comp.F, F), projective)
        for s, values, projective, _ in comp.rows
    ]
    found = structural_match(
        book_classes, book_rows, comp.invariants, comp_rows, F, None
    )
    assert found is not None, (
        f"no class and irrep bijection reproduces the server's double-group table of {symbol}"
    )
    _, pairing, _ = found
    assert sorted(pairing) == list(range(len(comp_rows)))
    for row, s in zip(data["irreps"], pairing, strict=True):
        assert REALITY_INDICATOR[row["reality"]] == comp.rows[s][3], (
            row["label"],
            comp.rows[s][0],
        )


@pytest.mark.reference
def test_reference_files_are_present() -> None:
    if not ALTMANN_HERZIG.is_dir() and not BILBAO.is_dir():
        pytest.skip(f"no transcriptions under {REFERENCE}")
    assert {"D3", "I", "Ih"} <= set(_files(ALTMANN_HERZIG))
    assert {"D3", "Oh"} <= set(_files(BILBAO))
    assert {"D3", "Oh"} <= set(_files(BILBAO, "-double"))


# ---- REF-3 and REF-5: the labelled comparison with Altmann & Herzig ----------------------------

AH_LABELS = REFERENCE / "altmann-herzig" / "labels"
AH_CONVENTION = "altmann-herzig"
# the character transcription spells these class labels differently from the
# labels transcription (whose README fixes ``sigma_x``) and from the reference
# class label of the convention (``15sigma`` of ``Ih``, docs/spec/conjugacy-classes.md)
BOOK_CLASS_SPELLING = {
    "sigmax": "sigma_x",
    "sigmay": "sigma_y",
    "sigmaz": "sigma_z",
    "15sigma_h": "15sigma",
}


def _labelled(symbol: str) -> dict[str, Any]:
    return model_document(symbol, label_convention=AH_CONVENTION)


@pytest.mark.reference
@pytest.mark.parametrize("symbol", _files(ALTMANN_HERZIG))
def test_ref5_reference_irrep_symbols_are_the_books(symbol: str) -> None:
    data = _read_altmann_herzig(symbol)
    ours = [
        record["reference_symbol"]
        for record in _labelled(symbol)["characters"]["irreps"]
    ]
    assert sorted(ours) == sorted(row["label"] for row in data["irreps"])


@pytest.mark.reference
@pytest.mark.parametrize(
    "path",
    sorted(AH_LABELS.glob("*.toml")) if AH_LABELS.is_dir() else [],
    ids=lambda path: path.stem,
)
def test_ref5_reference_irrep_symbols_match_the_symbol_row(path: Path) -> None:
    transcription = tomllib.loads(path.read_text())
    ours = [
        record["reference_symbol"]
        for record in _labelled(transcription["group"])["characters"]["irreps"]
    ]
    assert sorted(ours) == sorted(transcription["book"]["irreps"])


def _column_of(columns: list[dict[str, Any]], book_label: str) -> int:
    """The class of ``G-dagger`` the book's column names: the class whose reference
    label is the book's (the untilded half of a split class) or that merged
    class's unbarred half."""
    label = BOOK_CLASS_SPELLING.get(book_label, book_label)
    matches = [
        c
        for c, column in enumerate(columns)
        if column["reference_label"] == label
        or column["reference_label"].startswith(label + "+")
    ]
    assert len(matches) == 1, (book_label, [c["reference_label"] for c in columns])
    return matches[0]


@pytest.mark.reference
@pytest.mark.parametrize("symbol", _files(ALTMANN_HERZIG))
def test_ref3_labelled_comparison_of_the_character_tables(symbol: str) -> None:
    data = _read_altmann_herzig(symbol)
    document = _labelled(symbol)
    table = document["characters"]
    F_ours = cyclotomic_field(document["field"]["conductor"])
    F_book = cyclotomic_field(data["conductor"])
    F = cyclotomic_field(math.lcm(F_book.conductor, F_ours.conductor))
    book_rows = _book_rows(data, F_book, F)
    labels = [c["label"] for c in data["classes"]]
    for entry in data.get("known_discrepancy", []):
        book_rows = _swap_projective(book_rows, labels, *entry["swap_projective"])
    columns = [_column_of(table["classes"], label) for label in labels]
    row_of = {record["reference_symbol"]: k for k, record in enumerate(table["irreps"])}
    assert len(row_of) == len(table["irreps"]) == len(book_rows)
    mismatches = []
    for label, values, projective in book_rows:
        k = row_of[label]
        assert table["irreps"][k]["projective"] == projective, label
        ours = [
            F_ours.lift(F_ours.from_exact(table["exact"][k][c]), F) for c in columns
        ]
        mismatches += [
            (label, class_label, F.exact(mine), F.exact(theirs))
            for class_label, mine, theirs in zip(labels, ours, values, strict=True)
            if mine != theirs
        ]
    assert mismatches == [], "(irrep, class, ours, book)"
