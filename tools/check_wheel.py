"""Check a built wheel against the packaging acceptance criteria (wayfinder ticket #35).

    uv build && uv run python tools/check_wheel.py dist/pointgroup_tables-*.whl

The wheel must declare no runtime dependency and Python >= 3.13, ship the
convention data and the schema under ``pointgroup_tables/data/``, and ship
nothing of the prototypes, the tools, the tests or the reference
transcriptions (docs/spec/cli.md, "Size, dependencies and distribution").
"""

from __future__ import annotations

import sys
import zipfile
from pathlib import Path

REQUIRED = (
    "pointgroup_tables/__init__.py",
    "pointgroup_tables/data/schema.json",
    "pointgroup_tables/data/conventions/labels/altmann-herzig.toml",
    "pointgroup_tables/data/conventions/basis/canonical.toml",
    "pointgroup_tables/data/conventions/basis/altmann-herzig.toml",
    "pointgroup_tables/data/conventions/bethe/koster.toml",
)
FORBIDDEN_PREFIXES = ("prototypes/", "tools/", "tests/", "docs/")
FORBIDDEN_PARTS = ("reference", "samples", "__pycache__")


def check(wheel: Path) -> list[str]:
    problems: list[str] = []
    with zipfile.ZipFile(wheel) as archive:
        names = archive.namelist()
        metadata_name = next(n for n in names if n.endswith(".dist-info/METADATA"))
        metadata = archive.read(metadata_name).decode()
    for line in metadata.splitlines():
        if line.startswith("Requires-Dist:"):
            problems.append(f"runtime dependency declared: {line}")
    if "Requires-Python: >=3.13" not in metadata:
        problems.append("Requires-Python is not '>=3.13'")
    if "Name: pointgroup-tables" not in metadata:
        problems.append("the distribution is not named pointgroup-tables")
    for required in REQUIRED:
        if required not in names:
            problems.append(f"missing from the wheel: {required}")
    for name in names:
        if name.startswith(FORBIDDEN_PREFIXES) or any(
            part in FORBIDDEN_PARTS for part in name.split("/")
        ):
            problems.append(f"must not ship: {name}")
        if not (name.startswith("pointgroup_tables/") or ".dist-info/" in name):
            problems.append(f"outside the package: {name}")
        if name.endswith((".pyc", ".pyo")):
            problems.append(f"bytecode in the wheel: {name}")
    return problems


def main(argv: list[str]) -> int:
    if len(argv) != 1:
        print("usage: check_wheel.py <wheel>", file=sys.stderr)
        return 2
    wheel = Path(argv[0])
    problems = check(wheel)
    for problem in problems:
        print(problem)
    if problems:
        print(f"{wheel.name}: {len(problems)} problem(s)")
        return 1
    print(f"{wheel.name}: no runtime dependency, Python >= 3.13, package only")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
