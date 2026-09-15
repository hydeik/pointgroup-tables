# Altmann & Herzig symmetrized bases

Transcription of the tables T *n*.6 (*Symmetrized bases*) of S. L. Altmann
and P. Herzig, *Point-Group Theory Tables*, 2nd ed., for wayfinder ticket
[Transcribe Altmann & Herzig symmetrized bases and remaining axial matrices to verify the spin-harmonic rule](https://github.com/hydeik/pointgroup-tables/issues/22).
They are the book's own account of the bases behind its axial matrix tables
(T *n*.7) and verify the **spin-harmonic rule** of the `altmann-herzig`
basis convention (`docs/spec/basis-conventions.md`, ADR-0008) for the whole
range the book tabulates, `n <= 10`.

One TOML file per group, named by the normalised Schoenflies symbol; the
second setting of `C3v` is `C3vB.toml` (as in `../labels/`). 28 files:
`D5`, `D8`, `D9`, `D10` (T 25.6, 28.6, 29.6, 30.6), the masters `D2h` ...
`D10h` (T 31.6 - 39.6), `D2d`, `D4d`, `D6d`, `D8d`, `D10d` (T 41.6, 43.6,
45.6, 47.6, 49.6: the `Dnd` groups whose principal element is `S_2n`) and
`C2v` ... `C10v` (T 50.6 - 58.6, both settings of T 51.6).

```
uv run python tools/ah_axial_rule.py check    # every ket against the rule, every matrix table against both
uv run python tools/ah_axial_rule.py survey   # one line per ket: the book's (j, m1, s) next to the rule's (m1, m~)
```

## What is transcribed

Only the two-dimensional irreps (ordinary and projective; the book's
one-dimensional `1E`/`2E` pairs and their `|j m>+-` bases are not needed).
For every printed row of such an irrep the file records the **first ket**
of the row, `<|j m1>, s |j -m1>|`: `j`, `m1` (negative for the book's
overlined `m`), the sign `s` of the second partner, and whether the row is
marked with the book's bullet (`|j m>*`, an ungerade spin harmonic). The
row's further kets (higher `j`, the same `m1` up to the table's `+-mu`)
are not transcribed. The projective irreps of the improper groups print
two rows (gerade and ungerade harmonics), the odd `Cnv` groups four; every
row is kept, so a file lists more kets than irreps.

The OCR layer of the PDF drops the overline that marks a negative `m` and
often the minus sign of the second partner, so every table was read from
the page image (`get_pixmap`, zoom 2.2 for the page, 6 for doubtful rows).
Verify against the print.

## File layout

```toml
group = "D9h"
table = "T 38.6, p. 321"
status = "done"
source = "..."
transcribed_by = "..."
notes = ""

[irreps."E5/2"]                 # the book's irrep label
kets = [                        # one entry per printed row: the first ket <|j m1>, s |j -m1>|
  { j = "5/2", m1 = "-5/2", s = 1, bullet = false },
  { j = "13/2", m1 = "13/2", s = 1, bullet = true },
]
```

`setting = "A"` / `"B"` marks the two `C3v` tables.

## How the check reads a ket

`tools/ah_axial_rule.py` computes, exactly in `Q(zeta_M)`, the Wigner
matrices of the two harmonics under the group's generators: the principal
rotation `C_N` about `z` (`N = 2n` for `Dnh` with odd `n` and for `Dnd`,
whose principal element is `i C_2n`) and the rotation by `pi` about the pole
of the anchored secondary element (`x`, or `y` where the book's `sigma_v1`
of `Cnv` lies in the `yz` plane, i.e. `4` does not divide `n`). An
improper generator is the parity of the harmonic (`(-1)^j` for integral
`j`, `+1`/`-1` for a plain/bulleted spin harmonic) times the rotation's
matrix. The check then requires, for every ket:

- the exponent the ket carries on the principal element is the one the
  book's label names (`E_k`, `E_{k/2}`, the `'`/`''` of `Dnh` with odd `n`,
  the `S_2n`-based numbering of `Dnd`), and for `Dnh` with even `n` the
  harmonic's parity is the label's `g`/`u`;
- both generator matrices equal the rule's `D(a)`, `D(b)` entry by entry.

The generator matrices transcribed from T *n*.7 (`../matrices/`) are
compared with the rule by the same tool, so kets, printed matrices and
rule agree pairwise.

## Findings

Every ket of every table above agrees with the rule as amended by this
transcription (`docs/spec/basis-conventions.md`, "Spin-harmonic rule"):
the chain of levels divides `N` by its smallest prime factor at each step
(`3, 9` for `n = 9`, not `9` alone), the phase of the secondary element
comes from the level that decides the first partner, `E_{N/4}` of `D6d` and
`D10d` takes the negative partner first, and the ordinary irreps of the odd
`Cnv` groups carry the `D_n` matrices of `C2'(k)`. No printed ket had to be
corrected. One oddity is recorded in `C9v.toml`: the second and third rows
of `E5/2` are printed identically (`<|13/2 -5/2>, |13/2 5/2>|`), both kept
as printed; they are consistent with the first row.

## Redistribution

Committed to the repository as test data by the decision of
[Verification strategy](https://github.com/hydeik/pointgroup-tables/issues/12)
(ADR-0006). See `tests/reference/README.md` for the licensing note.
