"""``--verify``: the self-consistency checks of the document about to be emitted
(docs/spec/verification.md, "Where checks run").

The checkers read the tables the request produces plus the base tables they
need (the element list, the Cayley table and the factor table, and the
character table under ``matrices``), all built under the request's options
with every value form present, so that the exact values are there to verify
the numeric and pretty ones. The bytes emitted afterwards come from
``build_document`` and are the same with and without the flag.
"""

from __future__ import annotations

from dataclasses import replace
from typing import Any

from .checks import VERIFY_CHECKS, Failure, check
from .document import TABLE_ORDER, TABLES_OF_SUBCOMMAND, VALUE_FORMS, Request, _assemble

# the base tables every verification reads
BASE_TABLES = ("elements", "cayley", "factor")


def verification_tables(subcommand: str) -> tuple[str, ...]:
    """The tables of the subcommand's document plus the base tables, in table order."""
    wanted = set(TABLES_OF_SUBCOMMAND[subcommand]) | set(BASE_TABLES)
    if "matrices" in wanted:
        wanted.add("characters")  # MX-4 and IL-* read the characters
    return tuple(name for name in TABLE_ORDER if name in wanted)


def verification_checks(request: Request) -> list[str]:
    """The check identifiers ``--verify`` runs for a request, in ledger order:
    ``AG-*`` and ``CL-*`` always, ``DG-*`` when the double group is involved,
    ``CH-*`` with ``characters``, ``MX-*`` with ``matrices``, and ``IL-*`` and
    ``VL-*`` with either."""
    tables = TABLES_OF_SUBCOMMAND[request.subcommand]
    families = {"AG", "CL"}
    irreps = "characters" in tables or "matrices" in tables
    if request.double or (irreps and request.irreps != "ordinary"):
        families.add("DG")
    if "characters" in tables:
        families.add("CH")
    if "matrices" in tables:
        families.add("MX")
    if irreps:
        families.update(("IL", "VL"))
    return [name for name in VERIFY_CHECKS if name.split("-")[0] in families]


def verification_document(request: Request) -> dict[str, Any]:
    """The document the checkers read for a request."""
    model_request = replace(request, values=VALUE_FORMS)
    return _assemble(model_request, verification_tables(request.subcommand))


def verify(request: Request) -> list[Failure]:
    """Every failing check of the request's document, in ledger order."""
    return check(verification_document(request), verification_checks(request))
