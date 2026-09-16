# Koster Bethe-symbol transcription

Reference data for the `bethe` field of the irrep record
(`docs/spec/irrep-labels.md`, "Bethe aliases"), transcribed by hand from
G. F. Koster, J. O. Dimmock, R. G. Wheeler and H. Statz, *Properties of the
Thirty-Two Point Groups* (MIT Press, 1963), for wayfinder ticket
[Transcribe Koster Γ_i numbering for the 32 crystallographic groups](https://github.com/hydeik/pointgroup-tables/issues/15).
One TOML file per group, named by the normalised Schoenflies symbol
(`Cs` is `C1h.toml`, `Ci` is `S2.toml`). Templates are produced and checked by

```
uv run python tools/koster_bethe.py templates
uv run python tools/koster_bethe.py check            # every file
uv run python tools/koster_bethe.py check D3h Oh     # some groups
```

The shipped Koster table, `src/pointgroup_tables/data/conventions/bethe/koster.toml`,
is written from these files by

```
uv run python tools/koster_bethe.py generate
```

which reruns the checker, rebuilds every group with the package and refuses
to write when a transcribed symbol or character is not the one the package
computes. The files themselves are test data (`docs/spec/verification.md`,
`REF-8`), read by `tests/test_reference_bethe.py`.

## File layout

```toml
group = "D3h"
table = "Table 55, p. 71"     # where in the book
status = "done"               # todo | partial | done
notes = ""                    # disagreements with Bradley & Cracknell or Bilbao; anything unreadable
elements = ["S3^-1", "C3^1", "C2'(1)", "sigma_h"]   # identifying elements, canonical labels

[[irreps]]
symbol = "E''"                # canonical Mulliken symbol, canonical order; do not edit
dimension = 2
projective = false
characters = ["1", "-1", "0", "-2"]   # on `elements`; lifted elements for a projective irrep
bethe = "Gamma5"              # the one field to fill
```

Everything except `table`, `status`, `notes` and `bethe` is generated; the
checker rejects a file whose generated part no longer matches the computed
irreps (rerun `templates --force` after a rule change).

## How to fill a file

1. Open Koster's table for the group. Its rows are `Γ_1, Γ_2, ...`, ordinary
   irreps first, then the extra (double-group) irreps; its columns are the
   book's classes.
2. For each `[[irreps]]` entry, find the row of the same dimension whose
   characters on the listed elements equal `characters` (the `numeric`
   line, present when a value is not rational, gives the decimal form;
   `zeta3` is `exp(2πi/3) = ω`, `zeta6` is `exp(2πi/6)`, `i` is the
   imaginary unit, `sqrt2`, `sqrt3` are the positive roots). Only the
   identifying elements matter; the book's other columns can be ignored.
3. Write the row's symbol as `bethe`: `"Gamma5"`; in a centrosymmetric group
   `"Gamma5+"` or `"Gamma5-"` (the sign is Koster's parity superscript; it
   must agree with the `g`/`u` suffix of the canonical symbol).
4. A complex-conjugate pair (`1E`/`2E`, `1E1/2`/`2E1/2`, ...) is two rows in
   Koster too; keep the members apart by the sign of the imaginary part on
   the principal element (`1` has the positive exponent, i.e. `zeta3`, not
   `zeta3^2`).
5. Fill `table`, set `status`, and put into `notes` every place where
   Bradley & Cracknell (Table 2.2) or the Bilbao pages number the irrep
   differently; Koster decides.
6. Run the checker.

## Element mapping

The elements are canonical labels (`docs/spec/label-conventions.md`), not
Koster's. Koster orients the groups as Altmann & Herzig do in the cases
that matter here: `C2'(1)` is the two-fold axis along `x` (Koster's `C2'`
or `C2x`), `sigma_v(1)` the plane containing `z` and `x`, `C4^1(1)` the
four-fold axis along `z` in `O`/`Oh`, `S4^-1(1)` its improper counterpart
in `Td`, and `C3^1(1)` the three-fold axis along `[111]` rotating by
`+2π/3`. When a transcription depends on this orientation (the `1`/`2`
subscripts of `B` irreps in `D4`, `C4v`, `D2d`, the `B1`/`B2`/`B3` of
`D2`), say so in `notes` if the book's own axis statements were used to
decide.
