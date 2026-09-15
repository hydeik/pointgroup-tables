"""The group argument: Schoenflies grammar, aliases and the normalised symbol.

See docs/spec/cli.md ("The group argument") and docs/spec/abstract-groups.md
("Aliases").
"""

from __future__ import annotations

import re
from dataclasses import dataclass

from .errors import GROUPS_HINT, UnknownGroupError, UnsupportedOrderError, UsageError

GRAMMAR = re.compile(r"^([A-Za-z]+)(0|[1-9][0-9]*)?([A-Za-z]*)$")

AXIAL_FAMILIES = ("Cn", "Cnv", "Cnh", "S2n", "Dn", "Dnh", "Dnd")
POLYHEDRAL_FAMILIES = ("T", "Td", "Th", "O", "Oh", "I", "Ih")

# family key (the letters around n, canonical spelling) -> family name;
# axial keys take n, polyhedral keys forbid it
AXIAL_KEYS = {
    "C": "Cn",
    "Cv": "Cnv",
    "Ch": "Cnh",
    "S": "S2n",
    "D": "Dn",
    "Dh": "Dnh",
    "Dd": "Dnd",
}
_AXIAL_KEYS = {key.lower(): family for key, family in AXIAL_KEYS.items()}
_POLYHEDRAL_KEYS: dict[str, str] = {
    family.lower(): family for family in POLYHEDRAL_FAMILIES
}
# the alias table of docs/spec/abstract-groups.md, in its order
ALIASES = {
    "Cs": "C1h",
    "Ci": "S2",
    "C1v": "C1h",
    "D1": "C2",
    "D1h": "C2v",
    "D1d": "C2h",
    "V": "D2",
    "Vh": "D2h",
    "Vd": "D2d",
}
# the degenerate axial symbols among them (they carry an n) name another group
_AXIAL_ALIASES = {
    alias: symbol for alias, symbol in ALIASES.items() if alias[1:2].isdigit()
}
# the others (lower case) are matched as family keys without n
_ALIAS_KEYS = {
    alias.lower(): symbol
    for alias, symbol in ALIASES.items()
    if alias not in _AXIAL_ALIASES
}


@dataclass(frozen=True)
class GroupSpec:
    """A parsed group argument: the normalised symbol and its family."""

    symbol: str
    input_symbol: str
    family: str
    n: int | None

    @property
    def polyhedral(self) -> bool:
        return self.n is None

    @property
    def family_n(self) -> int:
        """The ``n`` of the family name: the subscript, or half of it for ``S2n``."""
        assert self.n is not None
        return self.n // 2 if self.family == "S2n" else self.n


def subscript(family: str, n: int) -> int:
    """The subscript of the symbol of an axial family at its ``n``: ``2n`` for ``S2n``."""
    return 2 * n if family == "S2n" else n


def axial_symbol(family: str, n: int) -> str:
    """The symbol of an axial family at axis order ``n`` (``S2n`` takes its subscript)."""
    head, tail = family[0], family[-1] if family[-1] in "vhd" else ""
    return f"{head}{n}{tail}"


def parse_group(argument: str) -> GroupSpec:
    """Parse a group argument as typed; raise the typed error of its exit code."""
    match = GRAMMAR.match(argument)
    if match is None:
        raise UsageError(
            f"group argument '{argument}' does not match the Schoenflies grammar",
            GROUPS_HINT,
        )
    before, digits, after = match.groups()
    key = (before + after).lower()
    if key in _ALIAS_KEYS and digits is None:
        return _resolved(_ALIAS_KEYS[key], argument)
    if key in _POLYHEDRAL_KEYS and digits is None:
        family = _POLYHEDRAL_KEYS[key]
        return GroupSpec(symbol=family, input_symbol=argument, family=family, n=None)
    if key in _AXIAL_KEYS and digits is not None:
        n = int(digits)
        family = _AXIAL_KEYS[key]
        if n < 1:
            raise UnsupportedOrderError(
                f"unknown group '{argument}': the axis order must be at least 1",
                GROUPS_HINT,
            )
        if family == "S2n" and n % 2:
            raise UnsupportedOrderError(
                f"unknown group '{argument}': S2n takes even subscripts only",
                GROUPS_HINT,
            )
        return _resolved(axial_symbol(family, n), argument)
    raise UnknownGroupError(f"unknown group '{argument}'", GROUPS_HINT)


def _resolved(symbol: str, input_symbol: str) -> GroupSpec:
    """Build the spec of a normalised symbol, resolving the degenerate axial aliases."""
    symbol = _AXIAL_ALIASES.get(symbol, symbol)
    if symbol in POLYHEDRAL_FAMILIES:
        return GroupSpec(
            symbol=symbol, input_symbol=input_symbol, family=symbol, n=None
        )
    match = GRAMMAR.match(symbol)
    assert match is not None and match.group(2) is not None, symbol
    before, digits, after = match.groups()
    family = _AXIAL_KEYS[(before + after).lower()]
    return GroupSpec(
        symbol=symbol, input_symbol=input_symbol, family=family, n=int(digits)
    )
