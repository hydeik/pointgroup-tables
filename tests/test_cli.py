"""The command line (docs/spec/cli.md, ``CLI-1``, ``CLI-2``, ``OS-6`` and the
utility subcommands)."""

from __future__ import annotations

import json
import tomllib
from pathlib import Path
from typing import Any

import pytest
import yaml
from test_symbols import GROUP_ARGUMENTS

from pointgroup_tables import main
from pointgroup_tables.formats import to_json
from pointgroup_tables.listing import shipped_schema


def run(
    capsysbinary: pytest.CaptureFixture[bytes], *argv: str
) -> tuple[int, bytes, str]:
    code = main(list(argv))
    captured = capsysbinary.readouterr()
    return code, captured.out, captured.err.decode()


@pytest.mark.parametrize(("argument", "code", "symbol"), GROUP_ARGUMENTS)
def test_cli1_group_argument_table(
    capsysbinary: pytest.CaptureFixture[bytes],
    argument: str,
    code: int,
    symbol: str | None,
) -> None:
    exit_code, out, err = run(capsysbinary, "elements", argument)
    assert exit_code == code
    if code == 0:
        document = json.loads(out)
        assert document["group"]["symbol"] == symbol
        assert document["group"]["input_symbol"] == argument
        assert err == ""
    else:
        assert out == b""
        assert "pointgroup-tables: error: " in err
        assert f"'{argument}'" in err
        if code == 2:
            assert err.startswith("usage: pointgroup-tables elements")
        else:
            assert err.endswith(
                "hint: run 'pointgroup-tables list groups' for the accepted symbols\n"
            )


@pytest.mark.parametrize(
    ("argv", "code"),
    [
        (["elements", "C3v", "--label-convention", "bilbao"], 4),
        (["elements", "C3v", "--basis-convention", "bilbao"], 4),
        (["factor", "C3v", "--double-group"], 2),
        (["elements", "C3v", "--irreps", "all"], 2),
        (["elements", "C3v", "--values", "float"], 2),
        (["elements", "C3v", "--format", "xml"], 2),
        (["elements"], 2),
        (["elements", "C3v", "D6h"], 2),
        (["classes", "C3v"], 2),
        ([], 2),
    ],
)
def test_cli2_failures_leave_the_stream_empty(
    capsysbinary: pytest.CaptureFixture[bytes],
    tmp_path: Path,
    argv: list[str],
    code: int,
) -> None:
    target = tmp_path / "out.json"
    exit_code, out, err = (
        run(capsysbinary, *argv, "-o", str(target)) if argv else run(capsysbinary)
    )
    assert exit_code == code
    assert out == b""
    assert "pointgroup-tables: error: " in err
    assert not target.exists()


def test_unknown_convention_hint(capsysbinary: pytest.CaptureFixture[bytes]) -> None:
    _, _, err = run(capsysbinary, "elements", "C3v", "--label-convention", "bilbao")
    assert err == (
        "pointgroup-tables: error: unknown label convention 'bilbao'\n"
        "hint: run 'pointgroup-tables list conventions' for the shipped conventions\n"
    )


def test_determinism(capsysbinary: pytest.CaptureFixture[bytes]) -> None:
    first = run(capsysbinary, "elements", "D6h", "--values", "exact")
    second = run(capsysbinary, "elements", "D6h", "--values", "exact")
    assert first == second
    assert first[0] == 0
    assert first[1].endswith(b"\n")
    assert b"\r" not in first[1]
    first[1].decode("ascii")


def test_output_file(
    capsysbinary: pytest.CaptureFixture[bytes], tmp_path: Path
) -> None:
    target = tmp_path / "table.yaml"
    target.write_text("stale")
    code, out, err = run(capsysbinary, "elements", "Oh", "-o", str(target))
    assert (code, out, err) == (0, b"", "")
    code, stdout, _ = run(capsysbinary, "elements", "Oh", "-o", "-")
    assert (
        target.read_bytes() == stdout
    )  # the format is never inferred from the file name
    assert json.loads(stdout)["tables"] == ["elements"]


def test_unwritable_output(
    capsysbinary: pytest.CaptureFixture[bytes], tmp_path: Path
) -> None:
    code, out, err = run(
        capsysbinary, "elements", "C1", "-o", str(tmp_path / "missing" / "out.json")
    )
    assert (code, out) == (1, b"")
    assert err.startswith("pointgroup-tables: error: cannot write ")


def test_version_and_help(capsysbinary: pytest.CaptureFixture[bytes]) -> None:
    code, out, _ = run(capsysbinary, "--version")
    assert code == 0
    assert out.decode().startswith("pointgroup-tables 0.") and out.count(b"\n") == 1
    code, out, _ = run(capsysbinary, "--help")
    assert code == 0
    assert b"usage: pointgroup-tables" in out
    code, out, _ = run(capsysbinary, "elements", "--help")
    assert code == 0
    assert b"--basis-convention" in out


def test_json_layout(capsysbinary: pytest.CaptureFixture[bytes]) -> None:
    code, out, _ = run(capsysbinary, "elements", "C3v")
    text = out.decode()
    assert code == 0
    assert '  "tables": ["elements"],\n' in text
    assert '        "members": [3, 4, 5],\n' in text
    assert text.endswith("}\n")
    assert json.loads(text) == json.loads(to_json(json.loads(text)))
    assert all(len(line) <= 100 for line in text.splitlines())


# ---- output formats, output files and the utility subcommands (ticket #26) ---

PARSERS = {"json": json.loads, "yaml": yaml.safe_load, "toml": tomllib.loads}


def parsed(data: bytes, fmt: str) -> Any:
    return PARSERS[fmt](data.decode("ascii"))


@pytest.mark.parametrize("fmt", ["json", "yaml", "toml"])
def test_os6_every_format_is_deterministic_and_writes_the_same_bytes_to_a_file(
    capsysbinary: pytest.CaptureFixture[bytes], tmp_path: Path, fmt: str
) -> None:
    runs = [
        run(capsysbinary, "elements", "D3h", "--format", fmt),
        run(capsysbinary, "elements", "D3h", "-f", fmt),
        run(capsysbinary, "elements", "D3h", "-f", fmt, "-o", "-"),
    ]
    assert runs[0] == runs[1] == runs[2]
    code, out, err = runs[0]
    assert (code, err) == (0, "")
    assert out.isascii() and out.endswith(b"\n") and b"\r" not in out
    target = tmp_path / "table.json"
    target.write_bytes(b"stale")
    assert run(capsysbinary, "elements", "D3h", "-f", fmt, "-o", str(target)) == (
        0,
        b"",
        "",
    )
    assert target.read_bytes() == out


@pytest.mark.parametrize("fmt", ["yaml", "toml"])
def test_yaml_and_toml_are_the_json_object(
    capsysbinary: pytest.CaptureFixture[bytes], fmt: str
) -> None:
    _, reference, _ = run(capsysbinary, "elements", "d4d", "--values", "both,pretty")
    _, out, _ = run(
        capsysbinary, "elements", "d4d", "--values", "both,pretty", "-f", fmt
    )
    assert parsed(out, fmt) == json.loads(reference)


@pytest.mark.parametrize("fmt", ["json", "yaml", "toml"])
def test_list_groups(capsysbinary: pytest.CaptureFixture[bytes], fmt: str) -> None:
    code, out, err = run(capsysbinary, "list", "groups", "-f", fmt)
    assert (code, err) == (0, "")
    document = parsed(out, fmt)
    assert list(document) == ["families", "polyhedral", "aliases"]
    assert [entry["key"] for entry in document["families"]] == [
        "C", "Cv", "Ch", "S", "D", "Dh", "Dd",
    ]  # fmt: skip
    assert list(document["families"][0]) == [
        "key", "family", "constraint", "rules", "example",
    ]  # fmt: skip
    assert [entry["symbol"] for entry in document["polyhedral"]] == [
        "T", "Td", "Th", "O", "Oh", "I", "Ih",
    ]  # fmt: skip
    assert {"input": "Vd", "symbol": "D2d"} in document["aliases"]
    assert run(capsysbinary, "list", "groups", "-f", fmt) == (code, out, err)


@pytest.mark.parametrize("fmt", ["json", "yaml", "toml"])
def test_list_conventions(capsysbinary: pytest.CaptureFixture[bytes], fmt: str) -> None:
    code, out, err = run(capsysbinary, "list", "conventions", "--format", fmt)
    assert (code, err) == (0, "")
    document = parsed(out, fmt)
    assert list(document) == ["label_conventions", "basis_conventions"]
    assert {"name": "canonical"} in document["basis_conventions"]
    assert all("source" in entry for entry in document["label_conventions"])
    assert run(capsysbinary, "list", "conventions", "--format", fmt) == (
        code,
        out,
        err,
    )


def test_list_writes_to_a_file(
    capsysbinary: pytest.CaptureFixture[bytes], tmp_path: Path
) -> None:
    target = tmp_path / "groups.toml"
    assert run(capsysbinary, "list", "groups", "-f", "toml", "-o", str(target)) == (
        0,
        b"",
        "",
    )
    _, out, _ = run(capsysbinary, "list", "groups", "-f", "toml")
    assert target.read_bytes() == out


def test_schema_prints_the_shipped_schema(
    capsysbinary: pytest.CaptureFixture[bytes], tmp_path: Path
) -> None:
    code, out, err = run(capsysbinary, "schema")
    assert (code, err) == (0, "")
    assert out == shipped_schema()
    assert json.loads(out)["$id"] == "urn:pointgroup-tables:schema:1"
    target = tmp_path / "schema.json"
    assert run(capsysbinary, "schema", "-o", str(target)) == (0, b"", "")
    assert target.read_bytes() == out


@pytest.mark.parametrize(
    "argv",
    [
        ["schema", "--format", "json"],
        ["schema", "-f", "yaml"],
        ["schema", "--values", "exact"],
        ["schema", "C3v"],
        ["list"],
        ["list", "irreps"],
        ["list", "groups", "--values", "exact"],
        ["list", "groups", "--label-convention", "x"],
        ["list", "conventions", "C3v"],
        ["list", "conventions", "--verify"],
        ["schema", "--verify"],
    ],
)
def test_utility_subcommands_reject_other_options(
    capsysbinary: pytest.CaptureFixture[bytes], tmp_path: Path, argv: list[str]
) -> None:
    target = tmp_path / "out"
    code, out, err = run(capsysbinary, *argv, "-o", str(target))
    assert (code, out) == (2, b"")
    assert err.startswith("usage: pointgroup-tables")
    assert "pointgroup-tables: error: " in err
    assert not target.exists()


def test_utility_help(capsysbinary: pytest.CaptureFixture[bytes]) -> None:
    code, out, _ = run(capsysbinary, "list", "--help")
    assert code == 0 and b"groups" in out and b"conventions" in out
    code, out, _ = run(capsysbinary, "list", "groups", "--help")
    assert code == 0 and b"--format" in out and b"--values" not in out
    code, out, _ = run(capsysbinary, "schema", "--help")
    assert code == 0 and b"--output" in out and b"--format" not in out


# ---- cayley, factor and --double-group (ticket #27) --------------------------


@pytest.mark.parametrize(
    ("argv", "tables", "double", "group"),
    [
        (["cayley", "D3"], ["elements", "cayley"], False, "G"),
        (["cayley", "D3", "--double-group"], ["elements", "cayley"], True, "G_double"),
        (["factor", "D3"], ["elements", "factor"], False, "G"),
        (["elements", "D3", "--double-group"], ["elements"], True, "G_double"),
    ],
)
def test_cayley_factor_and_double_group_documents(
    capsysbinary: pytest.CaptureFixture[bytes],
    argv: list[str],
    tables: list[str],
    double: bool,
    group: str,
) -> None:
    code, out, err = run(capsysbinary, *argv)
    assert (code, err) == (0, "")
    document = json.loads(out)
    assert document["tables"] == tables
    assert document["double"] is double
    assert document["elements"]["group"] == group
    order = 12 if double else 6
    assert len(document["elements"]["items"]) == order
    if len(tables) == 2:
        table = document[tables[1]]["table"]
        expected_group = group if tables[1] == "cayley" else "G"
        assert document[tables[1]]["group"] == expected_group
        assert len(table) == (order if tables[1] == "cayley" else 6)


def test_double_group_help_and_usage_error(
    capsysbinary: pytest.CaptureFixture[bytes],
) -> None:
    code, out, _ = run(capsysbinary, "cayley", "--help")
    assert code == 0 and b"--double-group" in out
    code, out, _ = run(capsysbinary, "factor", "--help")
    assert code == 0 and b"--double-group" not in out
    code, out, err = run(capsysbinary, "factor", "D3", "--double-group")
    assert (code, out) == (2, b"")
    assert err.startswith("usage: pointgroup-tables")
    assert "pointgroup-tables: error: " in err and "'--double-group'" not in err
    assert "--double-group" in err


# ---- characters -------------------------------------------------------------------


def test_characters_document_and_irreps_option(
    capsysbinary: pytest.CaptureFixture[bytes],
) -> None:
    code, out, err = run(capsysbinary, "characters", "D3", "--irreps", "projective")
    assert (code, err) == (0, "")
    document = json.loads(out)
    assert document["tables"] == ["elements", "characters"]
    assert document["irreps"] == "projective"
    assert [r["index"] for r in document["characters"]["irreps"]] == [3, 4, 5]
    assert document["characters"]["class_group"] == "G_double"
    code, out, err = run(capsysbinary, "characters", "D3", "--irreps", "spinor")
    assert (code, out) == (2, b"")
    assert "--irreps" in err
    code, out, _ = run(capsysbinary, "characters", "--help")
    assert code == 0 and b"--irreps" in out and b"--double-group" not in out
    code, out, err = run(capsysbinary, "characters", "D3", "--double-group")
    assert (code, out) == (2, b"")


def test_characters_of_a_polyhedral_group_are_not_available_yet(
    capsysbinary: pytest.CaptureFixture[bytes],
) -> None:
    code, out, err = run(capsysbinary, "characters", "Oh")
    assert (code, out) == (1, b"")
    assert err.startswith(
        "pointgroup-tables: error: the irreps of 'Oh' are not available"
    )


@pytest.mark.parametrize("fmt", ["json", "yaml", "toml"])
def test_characters_round_trip_through_every_format(
    capsysbinary: pytest.CaptureFixture[bytes], fmt: str
) -> None:
    code, out, err = run(
        capsysbinary, "characters", "D3h", "--format", fmt, "--values", "both,pretty"
    )
    assert (code, err) == (0, "")
    reference = json.loads(
        to_json(
            json.loads(
                run(capsysbinary, "characters", "D3h", "--values", "both,pretty")[1]
            )
        )
    )
    loaded = {"json": json.loads, "yaml": yaml.safe_load, "toml": tomllib.loads}[fmt](
        out.decode()
    )
    assert loaded == reference
