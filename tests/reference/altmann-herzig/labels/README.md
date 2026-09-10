# Altmann & Herzig element-label transcription

Reference data for the `altmann-herzig` label convention
(`docs/spec/label-conventions.md`), transcribed from S. L. Altmann and
P. Herzig, *Point-Group Theory Tables*, 2nd ed., for wayfinder ticket
[Transcribe Altmann & Herzig element labels](https://github.com/hydeik/pointgroup-tables/issues/13).
One TOML file per group, named by the normalised Schoenflies symbol
(`C3v.toml`, `D6h.toml`, `Oh.toml`; the book's `Ci` and `Cs` are `S2.toml`
and `C1h.toml`; the second setting of `C3v` is `C3vB.toml`).

```
uv run python tools/ah_book.py check [OCR_DIR]     # rebuild the book's groups, compare with the print
uv run python tools/ah_book.py generate            # (re)write every file here
uv run python tools/ah_transcription.py check      # validate the files, derive maps and swaps
```

## How the files were produced

The book parametrises every proper operation by a rotation angle `phi` and
a pole `n` (tables T *n*.1, book section 3-2), writes every improper
operation as `i g`, and declares the quaternion `[[cos(phi/2), sin(phi/2) n]]`
of `g` to be the untilded double-group element (book equation 3.30). The
parameters of the seven master tables (T 31.1 D2h, T 33.1 D4h, T 35.1 D6h,
T 37.1 D8h, T 39.1 D10h, T 71.1 Oh, T 75.1 Ih) and of C3v's own T 51.1A/B
are transcribed by hand in `tools/ah_book.py`, together with the
subgroup-element tables T *m*.0 (which operations each group contains and
how it names them) and the irrep symbol row of every character table
T *n*.4. Everything else is computed from the parameters by quaternion
multiplication and checked against the print: the reconstructed factor
tables agree entry for entry with the printed T 31.3, T 33.3, T 35.3,
T 37.3, T 39.3, T 71.3 and T 74.3 (9424 entries).

## File layout

```toml
group = "D3"
table = "T 23, p. 196"        # the group's own tables in the book
status = "done"
source = "S. L. Altmann and P. Herzig, ..."
transcribed_by = "agent (Claude) from the OCR text of the PDF; verify against the print"
notes = """..."""

[book]
elements = ["E", "C3+", "C3-", "C2'1", "C2'2", "C2'3"]   # book order, book spelling
master = "D6h"                                          # whose T m.1 / T m.0 the book refers to
irreps = ["A1", "A2", "E", "E1/2", "1E3/2", "2E3/2"]    # the book's irrep symbol row

[book.parameters]                 # phi and pole n of the rotation part, as printed
"C3+" = { phi = "2pi/3", n = ["0", "0", "1"] }
"C2'2" = { phi = "pi", n = ["-1/2", "sqrt3/2", "0"] }
"sigma_h" = { phi = "pi", n = ["0", "0", "1"], proper = "C2" }   # improper: i times `proper`

[map]                             # book label -> canonical label (derived, see below)
"C2'2" = "C2'(3)"

[multiplication]                  # reconstructed; "~" marks a tilded product
rows = [["E", "C3+", ...], ...]
```

`[map]` is derived from anchors (`E`, the principal rotation, the first
secondary axis or plane, `i`; for the polyhedral groups the second
generator is found by search) by closing the book's multiplication table
against the canonical Cayley table; every entry is verified as a
homomorphism by the checker. The checker also solves the coboundary between
the book's tilde marks and the canonical section and prints the
`[barred.<group>] swapped = [...]` list the label convention needs
(`docs/spec/double-groups.md`).

## Transliteration of the book's labels

| Printed | Written | Note |
| ------- | ------- | ---- |
| `C₃⁺`, `C₃⁻` | `C3+`, `C3-` | subscript inline, superscript sign trailing |
| `C₈³⁺`, `S₁₀³⁻`, `C₅²⁺` | `C8^3+`, `S10^3-`, `C5^2+` | power after `^` |
| `C₂′₁`, `C₂″₂` | `C2'1`, `C2''2` | primes as `'`; the axis index follows |
| `C₂ₓ`, `C₄ᵧ⁺`, `S₄ᶻ⁻` | `C2x`, `C4y+`, `S4z-` | Cartesian subscripts inline |
| `C₃₁⁺`, `S₆₁⁻`, `C₃,₁₀⁺`, `S₁₀,₁³⁻` | `C3,1+`, `S6,1-`, `C3,10+`, `S10,1^3-` | polyhedral axis index after a comma |
| `C₂′ₐ`, `C₂ₐ` | `C2'a`, `C2a` | |
| `σ_h`, `σ_v1`, `σ_d2`, `σ_x`, `σ_a` | `sigma_h`, `sigma_v1`, `sigma_d2`, `sigma_x`, `sigma_a` | |
| `E`, `i` | `E`, `i` | |
| `Ẽ`, `C̃₃⁺` (double group) | `~E`, `~C3+` | leading `~`; only inside `multiplication.rows` |
| `E₁/₂,g`, `¹E₃/₂`, `F₃/₂`, `A₃/₂` (irreps) | `E1/2,g`, `1E3/2`, `F3/2`, `A3/2` | as printed |

Parameter tokens: `0`, `1`, `1/2`, `sqrt3/2`, `1/sqrt2`, `1/sqrt3`, `c8`
(`cos pi/8`), `s8`, `c5`, `s5`, `c10`, `s10`, `cos_sigma`, `sin_sigma`,
`cos_rho'`, `sin_rho'` with `sigma = arctan((1 + sqrt5)/2)` and
`rho' = arctan((3 + sqrt5)/2)`, each with an optional leading `-`; angles
`0`, `pi/5`, `pi/4`, `pi/3`, `2pi/5`, `pi/2`, `3pi/5`, `2pi/3`, `3pi/4`,
`4pi/5`, `pi`.

## Facts recorded for the label convention

- **Poles.** The book assigns poles by its Rule 1 (section 3-2): the poles
  of a class are images of one another under the group (the `C2'` poles of
  `D6h` are the `C6` rotations of `(1 0 0)`, so `C2'2` has pole
  `(-1/2, sqrt3/2, 0)` and `C2'3` pole `(-1/2, -sqrt3/2, 0)`). There is no
  positive-hemisphere rule; the tilde assignment therefore differs from
  the canonical section on the elements the checker lists as `swapped`.
- **Anchors.** The `[map]` of `C2v` and `C6v` uses the letter-preserving
  anchors `sigma_y = sigma_v(1)` and `sigma_v1 = sigma_v(1)`, not the
  geometric `sigma_x` / `sigma_d1`, so that the book's `B1` is the canonical
  `B1` (ADR-0007, `docs/spec/label-conventions.md`).
- **Settings.** `C3v` has two settings: A (`sigma_v1` with pole `(0 1 0)`,
  the `sigma_v` planes of `D6h`, used inside `D3h`) and B (`sigma_v1` with
  pole `(1 0 0)`, the `sigma_d` planes of `D6h`, used inside `D3d`).
- **Irrep letters.** Four-dimensional irreps are `F` (`F` of `I`, `F3/2` of
  `O`, `1F3/2`/`2F3/2` of `T`); the real one-dimensional spinor irreps are
  `A1/2` (`C1`, `Ci`) and `A3/2` (`C3`, `S6`); spinor `g`/`u` is written
  after a comma (`E1/2,g`); `Cs` uses `A'`/`A''`.

## Redistribution

The map's driver decided (ticket #13) that these filled files are committed
to the repository. See `tests/reference/README.md` for the licensing note.
