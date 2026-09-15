"""The command line (docs/spec/cli.md): ``main(argv) -> exit code``."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from typing import Any, NoReturn

from . import conventions
from .document import (
    DEFAULT_VALUES,
    Request,
    build_document,
    parse_values,
    tool_version,
)
from .errors import PointGroupTablesError, UsageError
from .formats import EMITTERS
from .symbols import parse_group

PROG = "pointgroup-tables"
TABLE_SUBCOMMANDS = {
    "elements": "the element list and conjugacy classes of a group",
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


def _add_common_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "-f",
        "--format",
        choices=sorted(EMITTERS),
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
    return parser


def _request(namespace: argparse.Namespace) -> Request:
    return Request(
        subcommand=namespace.subcommand,
        group=namespace.group.input_symbol,
        values=namespace.values,
        label_convention=namespace.label_convention,
        basis_convention=namespace.basis_convention,
    )


def _report(error: PointGroupTablesError) -> int:
    print(f"{PROG}: error: {error.message}", file=sys.stderr)
    if error.hint:
        print(f"hint: {error.hint}", file=sys.stderr)
    return error.exit_code


def _write(data: bytes, output: str) -> None:
    if output == "-":
        sys.stdout.buffer.write(data)
        sys.stdout.buffer.flush()
    else:
        with open(output, "wb") as handle:
            handle.write(data)


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
        document = build_document(_request(namespace))
    except PointGroupTablesError as error:
        return _report(error)
    data = EMITTERS[namespace.format](document).encode("ascii")
    try:
        _write(data, namespace.output)
    except OSError as exc:
        return _report(
            PointGroupTablesError(f"cannot write '{namespace.output}': {exc.strerror}")
        )
    return 0
