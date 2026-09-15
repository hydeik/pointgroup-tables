"""The shipped label and basis conventions (docs/spec/label-conventions.md,
docs/spec/basis-conventions.md).

Version 1 ships the canonical basis and, once their tickets land, the
Altmann & Herzig conventions; a name outside these tables is exit code 4.
"""

from __future__ import annotations

from typing import Any

from .errors import CONVENTIONS_HINT, UnknownConventionError

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
