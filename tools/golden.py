"""The golden documents (docs/spec/verification.md, ``GOLD-1``): check and regeneration.

    uv run python tools/golden.py check         # every golden file equals the tool's bytes
    uv run python tools/golden.py regenerate    # rewrite tests/golden/ from the current tool

The golden set is fixed by the Verification chapter ("Test organisation"): the
groups ``C3v``, ``D6h``, ``Oh`` and ``Ih``, every table subcommand, the three
formats, default options, at ``tests/golden/<group>/<subcommand>.<format>``.
The files are the regression and determinism oracle of the test suite
(``tests/test_golden.py``); each one is produced by the command line end to
end, exactly as ``pointgroup-tables <subcommand> <group> -f <format> -o <path>``.

A regeneration is reviewed as a diff, never accepted blindly: every document
embeds ``generator.version``, so a release bump changes exactly one line of
every file, and any other hunk is a change to the tables that needs a reason
in the release notes (``docs/release.md``).
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

from pointgroup_tables import main

ROOT = Path(__file__).resolve().parent.parent
GOLDEN_DIR = ROOT / "tests" / "golden"
GROUPS = ("C3v", "D6h", "Oh", "Ih")
SUBCOMMANDS = ("elements", "cayley", "factor", "characters", "matrices", "all")
FORMATS = ("json", "yaml", "toml")


def golden_files() -> list[tuple[str, str, str, Path]]:
    """Every ``(group, subcommand, format, path)`` of the golden set, in a fixed order."""
    return [
        (group, subcommand, fmt, GOLDEN_DIR / group / f"{subcommand}.{fmt}")
        for group in GROUPS
        for subcommand in SUBCOMMANDS
        for fmt in FORMATS
    ]


def emit(group: str, subcommand: str, fmt: str, path: Path) -> None:
    """Run the command line end to end, writing the document to ``path``."""
    code = main([subcommand, group, "--format", fmt, "--output", str(path)])
    if code != 0:
        raise SystemExit(
            f"golden: '{subcommand} {group} --format {fmt}' exited with code {code}"
        )


def stray_files() -> list[Path]:
    """Files under the golden directory that are not part of the set."""
    expected = {path for _, _, _, path in golden_files()}
    return sorted(
        path
        for path in GOLDEN_DIR.rglob("*")
        if path.is_file() and path not in expected
    )


def regenerate() -> int:
    for group, subcommand, fmt, path in golden_files():
        path.parent.mkdir(parents=True, exist_ok=True)
        emit(group, subcommand, fmt, path)
        print(f"wrote {path.relative_to(ROOT)} ({path.stat().st_size} bytes)")
    for path in stray_files():
        print(f"stray file, not part of the golden set: {path.relative_to(ROOT)}")
    print("review the result with: git diff --stat tests/golden")
    return 0


def check() -> int:
    failures: list[str] = []
    with tempfile.TemporaryDirectory() as directory:
        for group, subcommand, fmt, path in golden_files():
            fresh = Path(directory) / f"{group}_{subcommand}.{fmt}"
            emit(group, subcommand, fmt, fresh)
            if not path.exists():
                failures.append(f"missing: {path.relative_to(ROOT)}")
            elif fresh.read_bytes() != path.read_bytes():
                failures.append(f"differs: {path.relative_to(ROOT)}")
    failures.extend(f"stray: {path.relative_to(ROOT)}" for path in stray_files())
    for line in failures:
        print(line)
    total = len(golden_files())
    if failures:
        print(f"golden: {len(failures)} problem(s) over {total} files")
        return 1
    print(f"golden: all {total} files equal the tool's output")
    return 0


def usage() -> int:
    print(__doc__.split("\n\n")[1] if __doc__ else "usage: golden.py check|regenerate")
    return 2


if __name__ == "__main__":
    commands = {"check": check, "regenerate": regenerate}
    command = commands.get(sys.argv[1]) if len(sys.argv) == 2 else None
    raise SystemExit(command() if command else usage())
