"""The golden documents (docs/spec/verification.md, ``GOLD-1``): every table
subcommand of ``C3v``, ``D6h``, ``Oh`` and ``Ih`` in every format with default
options equals the checked-in file byte for byte. The files are written by
``tools/golden.py regenerate`` and reviewed as a diff."""

from __future__ import annotations

from pathlib import Path

import pytest

from pointgroup_tables import main

GOLDEN_DIR = Path(__file__).parent / "golden"
GROUPS = ("C3v", "D6h", "Oh", "Ih")
SUBCOMMANDS = ("elements", "cayley", "factor", "characters", "matrices", "all")
FORMATS = ("json", "yaml", "toml")
GOLDEN_SET = [
    (group, subcommand, fmt)
    for group in GROUPS
    for subcommand in SUBCOMMANDS
    for fmt in FORMATS
]


@pytest.mark.parametrize(
    ("group", "subcommand", "fmt"),
    GOLDEN_SET,
    ids=[f"{g}-{s}-{f}" for g, s, f in GOLDEN_SET],
)
def test_gold1_document_equals_the_checked_in_file(
    capsysbinary: pytest.CaptureFixture[bytes], group: str, subcommand: str, fmt: str
) -> None:
    code = main([subcommand, group, "--format", fmt])
    captured = capsysbinary.readouterr()
    assert code == 0, captured.err
    golden = (GOLDEN_DIR / group / f"{subcommand}.{fmt}").read_bytes()
    assert captured.out == golden, (
        f"{subcommand} {group} --format {fmt} differs from the golden file; "
        "review the change and regenerate with 'uv run python tools/golden.py regenerate'"
    )


def test_the_golden_directory_holds_exactly_the_ledger_s_set() -> None:
    on_disk = sorted(
        path.relative_to(GOLDEN_DIR).as_posix()
        for path in GOLDEN_DIR.rglob("*")
        if path.is_file()
    )
    assert on_disk == sorted(f"{g}/{s}.{f}" for g, s, f in GOLDEN_SET)


def test_golden_files_embed_the_installed_version() -> None:
    """A release bump regenerates the set: the version line is the one hunk a
    bump is allowed to change."""
    from importlib import metadata

    version = metadata.version("pointgroup-tables").encode()
    for group, subcommand, fmt in GOLDEN_SET:
        text = (GOLDEN_DIR / group / f"{subcommand}.{fmt}").read_bytes()
        assert version in text
