"""The shipped label and basis conventions (docs/spec/label-conventions.md,
docs/spec/basis-conventions.md).

Version 1 ships the canonical basis and, once their tickets land, the
Altmann & Herzig conventions; a name outside these tables is exit code 4.
The polyhedral generator data of a basis convention (the images of the
lifted generator pair of ``2T``, ``2O``, ``2I`` in every irrep of dimension
two or more) is read from the convention's TOML file on first use.
"""

from __future__ import annotations

import tomllib
from dataclasses import dataclass
from functools import cache
from importlib import resources
from typing import Any

from .errors import CONVENTIONS_HINT, UnknownConventionError
from .field import Matrix, cyclotomic_field
from .values import parse

CANONICAL = "canonical"

# name -> metadata block ({ name, source }); reference conventions add coverage per group
LABEL_CONVENTIONS: dict[str, dict[str, Any]] = {}
BASIS_CONVENTIONS: dict[str, dict[str, Any]] = {CANONICAL: {"name": CANONICAL}}


def label_convention(name: str) -> dict[str, Any]:
    if name not in LABEL_CONVENTIONS:
        raise UnknownConventionError(
            f"unknown label convention '{name}'", CONVENTIONS_HINT
        )
    return dict(LABEL_CONVENTIONS[name])


def basis_convention(name: str) -> dict[str, Any]:
    if name not in BASIS_CONVENTIONS:
        raise UnknownConventionError(
            f"unknown basis convention '{name}'", CONVENTIONS_HINT
        )
    return dict(BASIS_CONVENTIONS[name])


# ---- polyhedral generator data ------------------------------------------------------


@dataclass(frozen=True)
class PolyhedralData:
    """The ``[polyhedral.<cover>]`` block of a basis convention: the conductor of
    the matrices, the canonical labels of the lifted generator pair (elements of
    the proper group ``T``, ``O`` or ``I``, lifted by the canonical section) and,
    per canonical irrep symbol, the images ``D(q1), D(q2)`` as matrices over
    ``Q(zeta_conductor)``."""

    cover: str
    conductor: int
    generators: tuple[str, str]
    irreps: dict[str, tuple[Matrix, Matrix]]


def _basis_file(name: str) -> dict[str, Any]:
    path = resources.files(__package__) / "data" / "conventions" / "basis"
    return tomllib.loads((path / f"{name}.toml").read_text(encoding="utf-8"))


@cache
def polyhedral_data(name: str, cover: str) -> PolyhedralData:
    """The polyhedral generator data of a basis convention for ``2T``, ``2O`` or ``2I``."""
    basis_convention(name)
    block = _basis_file(name)["polyhedral"][cover]
    F = cyclotomic_field(block["conductor"])
    irreps: dict[str, tuple[Matrix, Matrix]] = {}
    for symbol, entry in block["irreps"].items():
        first, second = entry["matrices"]
        irreps[symbol] = (
            [[parse(F, text) for text in row] for row in first],
            [[parse(F, text) for text in row] for row in second],
        )
    generators = block["generators"]
    return PolyhedralData(cover, F.conductor, (generators[0], generators[1]), irreps)
