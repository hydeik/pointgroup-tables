"""Which checks ``--verify`` runs for a request (docs/spec/verification.md,
"Where checks run"), and that its document carries what those checks read."""

from __future__ import annotations

import pytest

from pointgroup_tables import checks
from pointgroup_tables.document import Request
from pointgroup_tables.verify import (
    verification_checks,
    verification_document,
    verification_tables,
    verify,
)


def families(names: list[str]) -> set[str]:
    return {name.split("-")[0] for name in names}


def test_verify_checks_exclude_the_test_only_rows() -> None:
    assert set(checks.VERIFY_CHECKS) == set(checks.CHECKS) - {"CH-6", "VL-3"}
    assert list(checks.VERIFY_CHECKS) == [
        name for name in checks.CHECKS if name in checks.VERIFY_CHECKS
    ]  # ledger order


@pytest.mark.parametrize(
    ("request_", "expected"),
    [
        (Request("elements", "C3v"), {"AG", "CL"}),
        (Request("cayley", "C3v"), {"AG", "CL"}),
        (Request("elements", "C3v", double=True), {"AG", "CL", "DG"}),
        (Request("factor", "C3v"), {"AG", "CL"}),
        (Request("characters", "C3v"), {"AG", "CL", "DG", "CH", "IL", "VL"}),
        (
            Request("characters", "C3v", irreps="ordinary"),
            {"AG", "CL", "CH", "IL", "VL"},
        ),
        (
            Request("characters", "C3v", irreps="projective"),
            {"AG", "CL", "DG", "CH", "IL", "VL"},
        ),
        (Request("matrices", "C3v"), {"AG", "CL", "DG", "MX", "IL", "VL"}),
        (Request("matrices", "C3v", irreps="ordinary"), {"AG", "CL", "MX", "IL", "VL"}),
        (Request("all", "C3v"), {"AG", "CL", "DG", "CH", "MX", "IL", "VL"}),
        (
            Request("all", "C3v", double=True, irreps="ordinary"),
            {"AG", "CL", "DG", "CH", "MX", "IL", "VL"},
        ),
    ],
)
def test_the_check_families_of_a_request(request_: Request, expected: set[str]) -> None:
    names = verification_checks(request_)
    assert families(names) == expected
    assert names == [name for name in checks.VERIFY_CHECKS if name in names]
    assert "CH-6" not in names and "VL-3" not in names


@pytest.mark.parametrize(
    ("subcommand", "expected"),
    [
        ("elements", ("elements", "cayley", "factor")),
        ("cayley", ("elements", "cayley", "factor")),
        ("factor", ("elements", "cayley", "factor")),
        ("characters", ("elements", "cayley", "factor", "characters")),
        ("matrices", ("elements", "cayley", "factor", "characters", "matrices")),
        ("all", ("elements", "cayley", "factor", "characters", "matrices")),
    ],
)
def test_the_verification_document_holds_the_base_tables(
    subcommand: str, expected: tuple[str, ...]
) -> None:
    assert verification_tables(subcommand) == expected


@pytest.mark.parametrize("subcommand", ["elements", "characters", "matrices", "all"])
def test_the_verification_document_carries_every_value_form(subcommand: str) -> None:
    request = Request(
        subcommand, "D2d", values=("numeric",), double=subcommand in ("elements", "all")
    )
    document = verification_document(request)
    assert document["values"] == ["exact", "numeric", "pretty"]
    assert document["double"] == request.double
    for name in verification_tables(subcommand):
        assert name in document
    for name in verification_checks(request):
        for table in checks.CHECKS[name][1]:
            assert table in document, (name, table)
    assert verify(request) == []


def test_the_verification_document_keeps_the_request_options() -> None:
    document = verification_document(Request("characters", "D3h", irreps="ordinary"))
    assert document["irreps"] == "ordinary"
    assert document["characters"]["class_group"] == "G"
    assert all(not r["projective"] for r in document["characters"]["irreps"])
