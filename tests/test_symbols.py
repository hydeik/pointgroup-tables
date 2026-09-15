"""The group argument grammar (docs/spec/cli.md, "The group argument")."""

from __future__ import annotations

import pytest

from pointgroup_tables.errors import PointGroupTablesError
from pointgroup_tables.symbols import GroupSpec, parse_group

# the verification table of the command line chapter, plus the aliases
GROUP_ARGUMENTS = [
    ("C3v", 0, "C3v"),
    ("c3V", 0, "C3v"),
    ("Cs", 0, "C1h"),
    ("vd", 0, "D2d"),
    ("D1h", 0, "C2v"),
    ("S2", 0, "S2"),
    ("Ih", 0, "Ih"),
    ("3v", 2, None),
    ("C03", 2, None),
    ("C3-v", 2, None),
    ("Cnd", 3, None),
    ("T5", 3, None),
    ("C", 3, None),
    ("Xh", 3, None),
    ("Cs3", 3, None),
    ("Dnv", 3, None),
    ("C0", 5, None),
    ("D0h", 5, None),
    ("S3", 5, None),
    ("Ci", 0, "S2"),
    ("C1v", 0, "C1h"),
    ("D1", 0, "C2"),
    ("D1d", 0, "C2h"),
    ("V", 0, "D2"),
    ("VH", 0, "D2h"),
    ("s8", 0, "S8"),
    ("d12H", 0, "D12h"),
    ("td", 0, "Td"),
    ("C37", 0, "C37"),
]


@pytest.mark.parametrize(("argument", "code", "symbol"), GROUP_ARGUMENTS)
def test_group_argument_table(argument: str, code: int, symbol: str | None) -> None:
    if code == 0:
        spec = parse_group(argument)
        assert spec.symbol == symbol
        assert spec.input_symbol == argument
    else:
        with pytest.raises(PointGroupTablesError) as info:
            parse_group(argument)
        assert info.value.exit_code == code
        assert f"'{argument}'" in info.value.message


@pytest.mark.parametrize(
    ("argument", "expected"),
    [
        ("C3v", GroupSpec("C3v", "C3v", "Cnv", 3)),
        ("S8", GroupSpec("S8", "S8", "S2n", 8)),
        ("Cs", GroupSpec("C1h", "Cs", "Cnh", 1)),
        ("D1h", GroupSpec("C2v", "D1h", "Cnv", 2)),
        ("oh", GroupSpec("Oh", "oh", "Oh", None)),
    ],
)
def test_family_and_n(argument: str, expected: GroupSpec) -> None:
    assert parse_group(argument) == expected
