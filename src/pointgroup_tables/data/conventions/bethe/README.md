# Bethe symbols

The Koster table: the Bethe symbol `Gamma_i` of every irrep of the 32
crystallographic point groups, emitted as the `bethe` field of the irrep
record whatever label convention is selected (`docs/spec/irrep-labels.md`,
"Bethe aliases"). It is not a convention a flag selects.

| File | Source | Produced by |
| ---- | ------ | ----------- |
| `koster.toml` | G. F. Koster, J. O. Dimmock, R. G. Wheeler and H. Statz, *Properties of the Thirty-Two Point Groups* (1963) | `uv run python tools/koster_bethe.py generate`, from the transcriptions under `tests/reference/koster/bethe/` |

Generated files are never edited by hand: change the transcriptions and
regenerate. Each entry names its irrep by the characters on the group's
identifying elements; `tests/test_reference_bethe.py` (`REF-8`) proves the
attachment against the transcriptions.
