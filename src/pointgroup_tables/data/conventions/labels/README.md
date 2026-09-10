# Label conventions

One TOML file per label convention, selected with `--label-convention <name>`
(`docs/spec/label-conventions.md`). The canonical labelling has no file.

| File | Source | Produced by |
| ---- | ------ | ----------- |
| `altmann-herzig.toml` | S. L. Altmann and P. Herzig, *Point-Group Theory Tables*, 2nd ed. | `uv run python tools/ah_convention.py generate`, from the book's parameter tables and the transcriptions under `tests/reference/altmann-herzig/labels/` |

Generated files are never edited by hand: change the tool or the
transcriptions and regenerate. `tools/ah_convention.py check` proves that the
rules and the data reproduce the transcribed labels of every group.
