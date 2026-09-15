"""The command line (docs/spec/cli.md, ``CLI-1``, ``CLI-2`` and determinism)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from test_symbols import GROUP_ARGUMENTS

from pointgroup_tables import main
from pointgroup_tables.formats import to_json


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
        (["elements", "C3v", "--double-group"], 2),
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
