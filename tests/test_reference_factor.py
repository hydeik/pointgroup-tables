"""``REF-4``: the factor system against the tilde marks of the Altmann & Herzig
Cayley tables (docs/spec/verification.md, docs/spec/double-groups.md).

The book's tilde marks (the ``[multiplication]`` table of a labels
transcription) define ``omega_book`` on the book's labels; through the
``[map]`` of the transcription they become a factor system on the canonical
elements. The book's section differs from the canonical one on the elements
the label convention lists as ``swapped``: with ``eps = -1`` on those and
``+1`` elsewhere, ``omega_book(g, h) = eps(g) eps(h) eps(gh) omega(g, h)`` on
every pair. The polyhedral swapped lists are the ``[barred.<group>]`` tables of
the shipped convention; the axial ones are what the convention derives from
its ``[secondary]`` direction table (``tools/ah_convention.py check``), listed
here until the label-convention ticket ships that derivation.
"""

from __future__ import annotations

import tomllib
from pathlib import Path

import pytest
from conftest import model

REFERENCE = Path(__file__).resolve().parent / "reference" / "altmann-herzig" / "labels"
CONVENTION = (
    Path(__file__).resolve().parent.parent
    / "src"
    / "pointgroup_tables"
    / "data"
    / "conventions"
    / "labels"
    / "altmann-herzig.toml"
)
GROUPS = ["C2v", "D3", "D4", "D3h", "T", "O", "Ih"]
# the axial swapped lists derived from the convention's [secondary] table; the
# D3 list is also stated in docs/spec/double-groups.md ("Section overrides")
AXIAL_SWAPPED = {
    "C2v": ["sigma_d(1)"],
    "D3": ["C2'(2)", "C2'(3)"],
    "D4": ["C2''(2)"],
    "D3h": ["C2'(2)", "C2'(3)", "sigma_v(3)"],
}


def swapped_elements(symbol: str) -> list[str]:
    """The canonical labels whose Altmann & Herzig lift is our barred element."""
    if symbol in AXIAL_SWAPPED:
        return AXIAL_SWAPPED[symbol]
    barred = tomllib.loads(CONVENTION.read_text()).get("barred", {})
    return barred.get(symbol, {}).get("swapped", [])


def book_factor_system(symbol: str) -> tuple[list[list[int]], list[str]]:
    """``omega_book`` on canonical indices, and the canonical labels in index order."""
    path = REFERENCE / f"{symbol}.toml"
    if not path.is_file():
        pytest.skip(f"missing reference file {path}")
    data = tomllib.loads(path.read_text())
    if "multiplication" not in data or "map" not in data:
        pytest.skip(f"{path.name} has no [multiplication] or [map] table")
    labels = [record["label"] for record in model(symbol)["elements"]["items"]]
    canonical = {book: labels.index(label) for book, label in data["map"].items()}
    elements = data["book"]["elements"]
    rows = data["multiplication"]["rows"]
    n = len(labels)
    assert len(elements) == n and len(rows) == n
    omega_book = [[1] * n for _ in range(n)]
    products = [[0] * n for _ in range(n)]
    for x, row in enumerate(rows):
        for y, entry in enumerate(row):
            g, h = canonical[elements[x]], canonical[elements[y]]
            omega_book[g][h] = -1 if entry.startswith("~") else 1
            products[g][h] = canonical[entry.lstrip("~")]
    cayley = model(symbol)["cayley"]["table"]
    assert products == cayley, "the [map] is not an isomorphism onto the Cayley table"
    return omega_book, labels


@pytest.mark.reference
@pytest.mark.parametrize("symbol", GROUPS)
def test_ref4_factor_system_matches_the_tilde_marks(symbol: str) -> None:
    omega_book, labels = book_factor_system(symbol)
    document = model(symbol)
    omega, cayley = document["factor"]["table"], document["cayley"]["table"]
    swapped = swapped_elements(symbol)
    assert set(swapped) <= set(labels), swapped
    eps = [-1 if label in swapped else 1 for label in labels]
    n = len(omega)
    mismatches = [
        (labels[g], labels[h])
        for g in range(n)
        for h in range(n)
        if omega_book[g][h] != eps[g] * eps[h] * eps[cayley[g][h]] * omega[g][h]
    ]
    assert mismatches == [], "omega differs from the tilde marks after the swaps"


@pytest.mark.reference
def test_ref4_marks_are_present() -> None:
    present = [symbol for symbol in GROUPS if (REFERENCE / f"{symbol}.toml").is_file()]
    if not present:
        pytest.skip(f"no transcription under {REFERENCE}")
    assert present == GROUPS
