# Altmann & Herzig element-label transcription

Reference data for the `altmann-herzig` label convention
(`docs/spec/label-conventions.md`), transcribed by hand from
S. L. Altmann and P. Herzig, *Point-Group Theory Tables*, 2nd ed., for
wayfinder ticket
[Transcribe Altmann & Herzig element labels](https://github.com/hydeik/pointgroup-tables/issues/13).
One TOML file per group, named by the normalised Schoenflies symbol
(`C3v.toml`, `D6h.toml`, `Oh.toml`; the book's `Ci` and `Cs` are `S2.toml` and `C1h.toml`). Templates are produced and checked by

```
uv run python tools/ah_transcription.py templates
uv run python tools/ah_transcription.py check            # every file
uv run python tools/ah_transcription.py check D3 D3h     # some groups
```

The book's *labels* and *orientations* are what is transcribed; the mapping
onto canonical labels is derived from them and from the multiplication
table, so the files can be filled without knowing the canonical rules.

## File layout

```toml
group = "D3"
table = "T 20.1, p. 176"      # where in the book
status = "done"               # todo | partial | done
notes = ""                    # anything that could not be read unambiguously

[book]
elements = ["E", "C3+", "C3-", "C2'1", "C2'2", "C2'3"]       # book order, book spelling
operations = { "C3+" = "about z", "C2'1" = "about x" }        # the book's statements only
irreps = ["A1", "A2", "E", "E1/2", "1E3/2", "2E3/2"]          # the book's irrep symbol row

[map]                          # book label -> canonical label, the certain ones
"E" = "E"
"C3+" = "C3^1"
"C2'1" = "C2'(1)"

[multiplication]               # optional; see below
rows = [
  ["E",     "C3+",   "C3-",   "C2'1",  "C2'2",  "C2'3"],
  ["C3+",   "C3-",   "~E",    "C2'3",  "C2'1",  "~C2'2"],
  ...
]
```

## Transliteration of the book's labels

Copy the label as printed, then apply these rules so the file is plain
ASCII:

| Printed | Written | Note |
| ------- | ------- | ---- |
| `C₃⁺`, `C₃⁻` | `C3+`, `C3-` | subscript inline, superscript sign trailing |
| `C₂′₁`, `C₂″₂` | `C2'1`, `C2''2` | primes as `'`; the axis index follows |
| `C₂ₓ`, `C₂ᵧ`, `C₂ᶻ` | `C2x`, `C2y`, `C2z` | Cartesian subscripts inline |
| `σ_h`, `σ_v1`, `σ_d2`, `σ_x` | `sigma_h`, `sigma_v1`, `sigma_d2`, `sigma_x` | |
| `S₆⁻`, `S₄⁺` | `S6-`, `S4+` | |
| `C₃₁⁺`, `C₂₁` (polyhedral, axis index) | `C31+`, `C21` | subscript digits run together; note the axis in `operations` |
| `i`, `E` | `i`, `E` | |
| `Ẽ`, `C̃₃⁺` (double group) | `~E`, `~C3+` | leading `~`; only inside `multiplication.rows` |

A label that cannot be read unambiguously goes into `notes` with its
position; do not guess.

## What to fill, in order of value

1. **`book.elements`** for every group: all `|G|` labels in the book's order
   (the single group only; the tilded elements are implied).
2. **`book.operations`**: every axis or plane statement the book makes for
   this group (`C2'1` about `x`, `sigma_d1` containing `z` and bisecting
   `x`, `y`, ...). Only what the book states; no interpretation.
3. **`book.irreps`**: the irrep symbol row, as printed (`1E3/2`, `F3/2`,
   `A3/2`, ...).
4. **`map`**: the entries the axis statements make certain: `E`, the
   principal rotation (`C_n^+` maps to `C{n}^1`), the first secondary axis
   (`C2'1` or `sigma_v1` maps to `C2'(1)` or `sigma_v(1)`), `i`, `sigma_h`.
   The checker derives the rest when a multiplication table is present.
5. **`multiplication.rows`**, for the priority groups `C2v`, `D3`, `D4`,
   `D3h`, `T`, `O` first: the book's multiplication table with its tilde
   marks (`~` prefix on the product). This fixes the map completely and
   lets the checker compare the book's tilde assignment with the canonical
   section (`docs/spec/double-groups.md`).
6. **Hemisphere statement**: once, in `notes` of `C2v.toml`: how the book
   defines the positive hemisphere for binary rotations (which axis is
   tried first).

## Redistribution

The templates are committed. Whether *filled* files may be committed is an
open point of the map (licensing of transcribed reference tables). Until it
is settled, keep filled files out of commits; when it is, either commit them
or add a `.gitignore` here and record the decision in the ticket.
