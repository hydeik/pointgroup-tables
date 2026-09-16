"""The command line (docs/spec/cli.md): ``main(argv) -> exit code``."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from typing import Any, NoReturn

from . import conventions
from .document import (
    DEFAULT_VALUES,
    DOUBLE_GROUP_SUBCOMMANDS,
    IRREPS_CHOICES,
    IRREPS_SUBCOMMANDS,
    Request,
    build_document,
    parse_values,
    tool_version,
)
from .errors import PointGroupTablesError, UsageError, VerificationError
from .formats import EMITTERS
from .listing import conventions_document, groups_document, shipped_schema
from .symbols import parse_group
from .verify import verify

PROG = "pointgroup-tables"
TABLE_SUBCOMMANDS = {
    "elements": "the element list and conjugacy classes of a group",
    "cayley": "the Cayley table of a group",
    "factor": "the factor table of the spin double group of a group",
    "characters": "the character table of the ordinary and projective irreps of a group",
    "matrices": "the irrep matrices of a group in the selected basis convention",
    "all": "the union of the five tables of a group under one metadata block",
}
LIST_SUBCOMMANDS = {
    "groups": "what the group argument accepts",
    "conventions": "the shipped label and basis conventions",
}


class Parser(argparse.ArgumentParser):
    """argparse with the tool's error line and exit code 2 on usage errors."""

    def error(self, message: str) -> NoReturn:
        self.print_usage(sys.stderr)
        print(f"{PROG}: error: {message}", file=sys.stderr)
        raise SystemExit(2)


def _group_argument(text: str) -> Any:
    """Parse the group as typed; grammar failures are usage errors of the parser."""
    try:
        return parse_group(text)
    except UsageError as exc:
        raise argparse.ArgumentTypeError(exc.message) from None


def _values_argument(text: str) -> tuple[str, ...]:
    try:
        return parse_values(text)
    except UsageError as exc:
        raise argparse.ArgumentTypeError(exc.message) from None


def _add_output_options(parser: argparse.ArgumentParser, formats: bool) -> None:
    """``-f`` (unless ``formats`` is off: ``schema`` is always JSON) and ``-o``."""
    if formats:
        parser.add_argument(
            "-f",
            "--format",
            choices=list(EMITTERS),
            default="json",
            help="output format (default: json)",
        )
    parser.add_argument(
        "-o",
        "--output",
        metavar="<path>",
        default="-",
        help="output file; '-' is standard output",
    )


def _add_common_options(parser: argparse.ArgumentParser) -> None:
    _add_output_options(parser, formats=True)
    parser.add_argument(
        "--values",
        metavar="<set>",
        type=_values_argument,
        default=DEFAULT_VALUES,
        help="comma-separated subset of exact, numeric, both, pretty (default: both)",
    )
    parser.add_argument(
        "--label-convention", metavar="<name>", help="a shipped label convention"
    )
    parser.add_argument(
        "--basis-convention",
        metavar="<name>",
        default=conventions.CANONICAL,
        help="canonical or a shipped basis convention (default: canonical)",
    )
    parser.add_argument(
        "--verify",
        action="store_true",
        help="run the self-consistency checks on the document before emitting it",
    )


def build_parser() -> Parser:
    parser = Parser(prog=PROG, description="Exact tables of the finite point groups.")
    parser.add_argument(
        "--version", action="version", version=f"{PROG} {tool_version()}"
    )
    subparsers = parser.add_subparsers(
        dest="subcommand", metavar="<subcommand>", required=True
    )
    for name, help_text in TABLE_SUBCOMMANDS.items():
        sub = subparsers.add_parser(
            name, help=help_text, description=f"Print {help_text}."
        )
        sub.add_argument(
            "group",
            metavar="<group>",
            type=_group_argument,
            help="a Schoenflies symbol",
        )
        _add_common_options(sub)
        if name in DOUBLE_GROUP_SUBCOMMANDS:
            sub.add_argument(
                "--double-group",
                action="store_true",
                help="emit the table of the spin double group instead",
            )
        if name in IRREPS_SUBCOMMANDS:
            sub.add_argument(
                "--irreps",
                choices=list(IRREPS_CHOICES),
                default="all",
                help="the irreps to include (default: all)",
            )
    lister = subparsers.add_parser(
        "list",
        help="list the accepted groups or the shipped conventions",
        description="List what the group argument or the convention options accept.",
    )
    listings = lister.add_subparsers(dest="listing", metavar="<what>", required=True)
    for name, help_text in LIST_SUBCOMMANDS.items():
        sub = listings.add_parser(
            name, help=help_text, description=f"List {help_text}."
        )
        _add_output_options(sub, formats=True)
    schema = subparsers.add_parser(
        "schema",
        help="the JSON Schema of the table documents",
        description="Print the JSON Schema of the table documents (always JSON).",
    )
    _add_output_options(schema, formats=False)
    return parser


def _request(namespace: argparse.Namespace) -> Request:
    return Request(
        subcommand=namespace.subcommand,
        group=namespace.group.input_symbol,
        values=namespace.values,
        label_convention=namespace.label_convention,
        basis_convention=namespace.basis_convention,
        double=getattr(namespace, "double_group", False),
        irreps=getattr(namespace, "irreps", "all"),
    )


def _report(error: PointGroupTablesError) -> int:
    print(f"{PROG}: error: {error.message}", file=sys.stderr)
    if error.hint:
        print(f"hint: {error.hint}", file=sys.stderr)
    for line in error.details:
        print(line, file=sys.stderr)
    return error.exit_code


def _write(data: bytes, output: str) -> None:
    if output == "-":
        sys.stdout.buffer.write(data)
        sys.stdout.buffer.flush()
    else:
        with open(output, "wb") as handle:
            handle.write(data)


def _render(namespace: argparse.Namespace) -> bytes:
    """The bytes of the invocation's document: nothing is written before this returns."""
    if namespace.subcommand == "schema":
        return shipped_schema()
    if namespace.subcommand == "list":
        document = (
            groups_document()
            if namespace.listing == "groups"
            else conventions_document()
        )
    else:
        request = _request(namespace)
        if namespace.verify:
            failures = verify(request)
            if failures:
                raise VerificationError(
                    namespace.group.symbol, tuple(str(f) for f in failures)
                )
        document = build_document(request)
    return EMITTERS[namespace.format](document).encode("ascii")


def main(argv: Sequence[str] | None = None) -> int:
    """Run the tool on ``argv`` (default: ``sys.argv[1:]``) and return the exit code."""
    parser = build_parser()
    try:
        namespace = parser.parse_args(argv)
    except SystemExit as exc:
        return exc.code if isinstance(exc.code, int) else 1
    except PointGroupTablesError as error:  # exit codes 3 and 5 from the group argument
        return _report(error)
    try:
        data = _render(namespace)
    except PointGroupTablesError as error:
        return _report(error)
    try:
        _write(data, namespace.output)
    except OSError as exc:
        return _report(
            PointGroupTablesError(f"cannot write '{namespace.output}': {exc.strerror}")
        )
    return 0
