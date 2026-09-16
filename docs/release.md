# Release checklist

A release is a version of the tool and of the conventions it ships: a
convention has no version of its own, and a document is reproducible from
`generator.version` alone (`docs/spec/cli.md`, "Size, dependencies and
distribution"). Every step below is run from a clean checkout of `main`.

## Before tagging

1. **The full suite is green**, every marker included, on Python 3.13:

   ```
   uv sync
   uv run ruff check . && uv run ruff format --check . && uv run ty check
   uv run pytest -m "" -rs
   ```

   `reference` tests that skip do so with a reason naming a missing file
   under `tests/reference/`; a skip for any other reason is a bug.
2. **Bump the version** in `pyproject.toml` (`[project] version`) and run
   `uv lock` so `uv.lock` records it. The version is read from the installed
   package metadata by `--version` and by every document's
   `generator.version`; nothing else needs editing.
3. **Regenerate the golden documents** and review the diff:

   ```
   uv sync
   uv run python tools/golden.py regenerate
   git diff --stat tests/golden
   ```

   A version bump changes exactly one line per file (the `generator.version`
   value). Any other hunk is a change to the tables: it must be intended,
   named in the release notes, and explained by a change to the code or the
   convention data since the last release. Never regenerate to make a failing
   `GOLD-1` pass without reading the diff.
4. **Convention data.** If a shipped convention file changed, the generating
   tool's `check` command still reproduces every transcription:

   ```
   uv run python tools/ah_convention.py check
   uv run python tools/ah_basis.py check
   uv run python tools/koster_bethe.py check
   ```

   and the change is noted in the release notes (a corrected transcription,
   a new rule, an extended `[coverage]` range), because a reader of an older
   document sees the difference only through `generator.version`.
5. **Build and inspect the wheel**:

   ```
   rm -rf dist && uv build
   uv run python tools/check_wheel.py dist/pointgroup_tables-*.whl
   ```

   The check requires: no `Requires-Dist`, `Requires-Python: >=3.13`, the
   convention data and `schema.json` present, nothing from `prototypes/`,
   `tools/`, `tests/` or the reference transcriptions.
6. **Install the wheel as a tool** and smoke-test it outside the checkout:

   ```
   uv tool install --reinstall dist/pointgroup_tables-*.whl
   pointgroup-tables --version
   pointgroup-tables characters C3v --verify
   pointgroup-tables matrices D6h --basis-convention altmann-herzig --verify > /dev/null
   pointgroup-tables list conventions
   uv tool uninstall pointgroup-tables
   ```

7. **Release notes.** One entry per user-visible change, grouped as: new
   groups or tables, convention data changes (with the affected groups and
   `coverage`), document changes (any key added; keys are never renamed or
   reordered within `schema_version = 1`), fixes. Reference misprints found
   and recorded as known discrepancies are listed too.
8. **Commit** the version bump, the lock file, the golden documents and the
   notes in one commit, and tag it `v<version>`.

## Publishing

9. `uv publish` (a PyPI API token in `UV_PUBLISH_TOKEN`), or upload the
   `dist/` artifacts that the `package` job of the CI workflow attached to
   the run of the tagged commit. Publish the sdist and the wheel together.
10. **Verify from PyPI** on a machine without the checkout:

    ```
    uv tool install pointgroup-tables==<version>
    pointgroup-tables --version
    pointgroup-tables all Oh --verify > /dev/null
    ```

11. Create the GitHub release from the tag with the release notes.

## Before the first public release

- Choose and add a licence file at the repository root and the matching
  `license` field and classifier in `pyproject.toml`. The reference
  transcriptions under `tests/reference/` are redistributed under the
  repository licence (ADR-0006), so the choice covers them too; keep the
  Bilbao citations of `tests/reference/README.md` in any redistribution.
- Confirm the author e-mail and the project URLs in `pyproject.toml`.
- Register the PyPI project name and a trusted publisher or an API token.
