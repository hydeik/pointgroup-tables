"""``REF-7``: under ``--basis-convention altmann-herzig`` the matrices of the
transcribed generators equal the book's entry by entry, and the Wigner
matrices of every ket of the symmetrized-bases transcription equal the
rule's (docs/spec/basis-conventions.md, docs/spec/verification.md).

Three transcriptions serve: the polyhedral matrix tables
(``tests/reference/altmann-herzig/matrices/{T,O,I}.toml``; every printed
operation of ``T`` and ``O``, the three generators of ``I`` with the two
corrected misprints of T 74.7b applied by the transcription), the ten axial
matrix tables (the same directory) and the 28 symmetrized-bases tables
(``.../bases/``). Book labels are matched to elements and irreps through the
reference labels of the ``altmann-herzig`` label convention (``REF-5``,
``REF-6``); a book operation whose untilded lift is our barred element
(``swapped``) carries ``-D(g)`` in the projective irreps.
"""

from __future__ import annotations

import functools
import tomllib
from fractions import Fraction
from pathlib import Path
from typing import Any

import pytest
from conftest import model as canonical_model

from pointgroup_tables import basis
from pointgroup_tables.abstract import point_group
from pointgroup_tables.document import model_document
from pointgroup_tables.field import CyclotomicField, Matrix, cyclotomic_field
from pointgroup_tables.irreps import named_elements
from pointgroup_tables.symbols import parse_group
from pointgroup_tables.values import parse

AH = "altmann-herzig"
REFERENCE = Path(__file__).resolve().parent / "reference" / "altmann-herzig"
MATRICES, BASES = REFERENCE / "matrices", REFERENCE / "bases"
POLYHEDRAL = ["T", "O", "I"]
AXIAL = ["D3", "D5", "D8", "D9", "D10", "C4v", "C6v", "D3h", "D2d", "D6h"]
BASES_FILES = sorted(BASES.glob("*.toml")) if BASES.is_dir() else []


def _load(path: Path) -> dict[str, Any]:
    if not path.is_file():
        pytest.skip(f"missing reference file {path}")
    return tomllib.loads(path.read_text())


@functools.cache
def document(symbol: str) -> dict[str, Any]:
    """The model under both conventions, in the index space of ``G-dagger`` so
    that the reference section (the tilde of every lift) can be read off."""
    return model_document(symbol, double=True, label_convention=AH, basis_convention=AH)


def _elements(doc: dict[str, Any]) -> dict[str, tuple[int, bool]]:
    """book label -> (index in G, whether the book's untilded operation is our barred lift)."""
    N = doc["group"]["order"]
    out = {}
    for record in doc["elements"]["items"][:N]:
        label = record["reference_label"]
        swapped = label.startswith("~")
        out[label[1:] if swapped else label] = (record["index"], swapped)
    return out


def _matrices(F: CyclotomicField, entry: dict[str, Any]) -> list[Matrix]:
    return [[[F.from_exact(x) for x in row] for row in D] for D in entry["exact"]]


def _book_matrix(
    F_book: CyclotomicField, F: CyclotomicField, rows: list[list[str]]
) -> Matrix:
    return [[F_book.lift(parse(F_book, text), F) for text in row] for row in rows]


def _ours(
    F: CyclotomicField, entry: dict[str, Any], index: int, swapped: bool
) -> Matrix:
    D = _matrices(F, entry)[index]
    if swapped and entry["projective"]:
        return [[F.neg(x) for x in row] for row in D]
    return D


# ---- the polyhedral tables ---------------------------------------------------------------


@pytest.mark.reference
@pytest.mark.parametrize("symbol", POLYHEDRAL)
def test_ref7_polyhedral_matrices_equal_the_book(symbol: str) -> None:
    data = _load(MATRICES / f"{symbol}.toml")
    doc = document(symbol)
    F, F_book = (
        cyclotomic_field(doc["field"]["conductor"]),
        cyclotomic_field(data["conductor"]),
    )
    assert doc["basis_convention"]["coverage"] == "reference"
    elements = _elements(doc)
    entries = {e["symbol"]: e for e in doc["matrices"]["irreps"]}
    canonical_of = data["canonical"]
    mismatches = []
    for book_symbol, block in data["irreps"].items():
        entry = entries[canonical_of[book_symbol]]
        for book_label, rows in zip(
            data["book"]["elements"], block["matrices"], strict=True
        ):
            index, swapped = elements[book_label]
            if _ours(F, entry, index, swapped) != _book_matrix(F_book, F, rows):
                mismatches.append((book_symbol, book_label))
    assert mismatches == []
    if symbol == "I":
        assert len(data["known_discrepancy"]) == 4  # the corrected misprints of T 74.7b
    else:
        assert "known_discrepancy" not in data


# ---- the axial tables --------------------------------------------------------------------


@pytest.mark.reference
@pytest.mark.parametrize("symbol", AXIAL)
def test_ref7_axial_generator_matrices_equal_the_book(symbol: str) -> None:
    data = _load(MATRICES / f"{symbol}.toml")
    doc = document(symbol)
    F, F_book = (
        cyclotomic_field(doc["field"]["conductor"]),
        cyclotomic_field(data["conductor"]),
    )
    assert F.conductor == F_book.conductor
    elements = _elements(doc)
    entries = {e["reference_symbol"]: e for e in doc["matrices"]["irreps"]}
    mismatches = []
    for book_symbol, block in data["irreps"].items():
        entry = entries[book_symbol]
        for book_label, rows in zip(data["elements"], block["matrices"], strict=True):
            index, swapped = elements[book_label]
            if _ours(F, entry, index, swapped) != _book_matrix(F_book, F, rows):
                mismatches.append((book_symbol, book_label))
    assert mismatches == []


# ---- the symmetrized bases --------------------------------------------------------------


def _kind(spec: Any) -> str:
    """``Dn``, ``Cnv``, ``Dnh_even``, ``Dnh_odd`` or ``Dnd_even``: where the parity of
    the harmonic attaches (to the secondary element of ``Cnv``, to the principal
    element ``i C_2n`` of ``Dnh`` with odd ``n`` and of ``Dnd``, to ``i`` of ``Dnh``
    with even ``n``)."""
    family, n = spec.family, spec.family_n
    if family == "Dnh":
        return "Dnh_even" if n % 2 == 0 else "Dnh_odd"
    if family == "Dnd":
        assert n % 2 == 0, spec.symbol
        return "Dnd_even"
    assert family in ("Dn", "Cnv"), spec.symbol
    return family


def _ket(
    F: CyclotomicField, N: int, j2: int, m2: int, s: int, bullet: bool, twisted: bool
) -> tuple[Matrix, Matrix, int]:
    """The Wigner matrices of the rotations ``C_N`` about ``z`` and ``C_2`` about the
    anchored pole (``x``, or ``y`` when ``twisted``) in the basis
    ``<|j m1>, s |j -m1>|`` (``j2 = 2j``, ``m2 = 2 m1``; Condon-Shortley harmonics,
    ``g f(r) = f(g^-1 r)``, ``exp(-i phi n.sigma/2)`` for the spin harmonics), and
    the parity of the harmonic (``(-1)^j`` for integral ``j``; ``-1`` for a bulleted
    spin harmonic)."""
    M = F.conductor
    assert M % (2 * N) == 0 and M % 4 == 0
    zero, sign = F.zero, F.rational(s)
    Da = [
        [F.zeta(-m2 * M // (2 * N) % M), zero],
        [zero, F.zeta(m2 * M // (2 * N) % M)],
    ]
    # C2 about the pole at azimuth alpha: |j m> -> exp(-i pi j) exp(2 i m alpha) |j -m>
    phase = F.zeta(-j2 * M // 4 % M)
    y = F.zeta(-m2 * M // 4 % M) if twisted else F.one
    Db = [
        [zero, F.mul(F.mul(phase, y), sign)],
        [F.mul(F.mul(phase, F.conjugate(y)), sign), zero],
    ]
    if j2 % 2 == 0:
        parity = -1 if (j2 // 2) % 2 else 1
    else:
        parity = -1 if bullet else 1
    return Da, Db, parity


def _scaled(F: CyclotomicField, D: Matrix, sign: int) -> Matrix:
    return D if sign == 1 else [[F.neg(x) for x in row] for row in D]


def _canonical_symbol(kind: str, n: int, N: int, label: str) -> str:
    """The canonical symbol of the irrep behind the book's label in a T n.6 table:
    the label itself (``E1/2,g`` spelled ``E1/2g``), except that ``E_k`` of ``Dnd``
    is numbered by ``C_n = S_2n^2`` in the book and by ``S_2n`` here (``k' = k (n + 1)``
    modulo ``2n``, reduced to ``(0, n)``)."""
    if kind == "Dnd_even" and "/" not in label:
        k = int(label[1:] or "1")
        kk = (k * (n + 1)) % (2 * n)
        kk = min(kk, 2 * n - kk)
        return "E" if N <= 4 else f"E{kk}"
    return label.replace(",", "")


@pytest.mark.reference
@pytest.mark.parametrize("path", BASES_FILES, ids=lambda path: path.stem)
def test_ref7_kets_of_the_symmetrized_bases_equal_the_rule(path: Path) -> None:
    """Every printed ket of every two-dimensional irrep gives, under the group's
    generators, the matrices the document carries for the irrep the book's label
    names (its principal element, its secondary element and, for ``Dnh`` with
    even ``n``, the inversion); the label is read as the tool that verified the
    rule reads it (``tools/ah_axial_rule.py``), by the exponent it names, not
    through the label convention's irrep overrides. The second setting of ``C3v`` (``sigma_v1`` on
    ``x``) is compared with the rule at ``alpha = 0``, which the shipped
    convention does not realise."""
    data = _load(path)
    symbol, setting = data["group"], data.get("setting")
    spec = parse_group(symbol)
    kind = _kind(spec)
    G = point_group(spec)
    named = named_elements(G)
    assert named.secondary is not None
    N = G.rotation_group.axis_order
    assert N is not None
    twisted = kind == "Cnv" and setting != "B" and spec.family_n % 4 != 0
    doc = document(symbol)
    F = cyclotomic_field(doc["field"]["conductor"])
    entries = {e["symbol"]: e for e in doc["matrices"]["irreps"]}
    problems = []
    count = 0
    for label, block in data["irreps"].items():
        entry = entries[_canonical_symbol(kind, spec.family_n, N, label)]
        ours = _matrices(F, entry)
        for row in block["kets"]:
            count += 1
            j2, m2 = int(2 * Fraction(row["j"])), int(2 * Fraction(row["m1"]))
            Ka, Kb, parity = _ket(F, N, j2, m2, row["s"], row["bullet"], twisted)
            where = f"{label} <|{row['j']} {row['m1']}>, {row['s']:+d}|{row['j']} -m1>|"
            if kind == "Cnv":
                Kb = _scaled(F, Kb, parity)
            elif kind in ("Dnh_odd", "Dnd_even"):
                Ka = _scaled(F, Ka, parity)
            elif kind == "Dnh_even":
                assert named.inversion is not None
                expected_parity = 1 if label.endswith("g") else -1
                if parity != expected_parity:
                    problems.append(
                        f"{where}: the parity of the harmonic is not the label's"
                    )
                if ours[named.inversion] != _scaled(
                    F, [[F.one, F.zero], [F.zero, F.one]], parity
                ):
                    problems.append(f"{where}: D(i) is not the parity of the harmonic")
            if setting == "B":
                k2 = 2 if label == "E" else 1
                Da, Db = basis.spin_harmonic_images(F, N, k2, 0)
            else:
                Da, Db = ours[named.principal], ours[named.secondary]
            if Da != Ka:
                problems.append(f"{where}: D(a) differs")
            elif Db != Kb:
                problems.append(f"{where}: D(b) differs")
    assert count > 0
    assert problems == []


# ---- presence ------------------------------------------------------------------------------


@pytest.mark.reference
def test_ref7_transcriptions_are_present() -> None:
    present = [
        symbol
        for symbol in POLYHEDRAL + AXIAL
        if (MATRICES / f"{symbol}.toml").is_file()
    ]
    if not present and not BASES_FILES:
        pytest.skip(f"no transcription under {REFERENCE}")
    assert present == POLYHEDRAL + AXIAL
    assert len(BASES_FILES) == 28
    assert {path.stem for path in BASES_FILES} >= {"C3v", "C3vB", "D6d", "D10d", "D9h"}
    assert canonical_model("D6h")["basis_convention"] == {"name": "canonical"}
