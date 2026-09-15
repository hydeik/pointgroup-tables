"""Typed errors that map to the exit codes of the command line chapter."""

from __future__ import annotations


class PointGroupTablesError(Exception):
    """Base of every error the tool reports as ``pointgroup-tables: error:``."""

    exit_code = 1

    def __init__(self, message: str, hint: str | None = None) -> None:
        super().__init__(message)
        self.message = message
        self.hint = hint


class UsageError(PointGroupTablesError):
    """Exit code 2: an argument that does not match what the tool accepts."""

    exit_code = 2


class UnknownGroupError(PointGroupTablesError):
    """Exit code 3: the argument matches the grammar but names no group."""

    exit_code = 3


class UnknownConventionError(PointGroupTablesError):
    """Exit code 4: a convention that is not shipped or does not cover the group."""

    exit_code = 4


class UnsupportedOrderError(PointGroupTablesError):
    """Exit code 5: a known family with an axis order out of range."""

    exit_code = 5


GROUPS_HINT = "run 'pointgroup-tables list groups' for the accepted symbols"
CONVENTIONS_HINT = (
    "run 'pointgroup-tables list conventions' for the shipped conventions"
)
