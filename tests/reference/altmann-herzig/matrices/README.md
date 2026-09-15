# Altmann & Herzig matrix representations

Transcriptions of the tables T *n*.7 of S. L. Altmann and P. Herzig,
*Point-Group Theory Tables*, 2nd ed. Two wayfinder tickets share this
directory: the axial samples of
[Transcribe Altmann & Herzig character tables and sample axial irrep matrices](https://github.com/hydeik/pointgroup-tables/issues/19)
(`D3`, `C4v`, `D3h`, `D2d`, `D6h`, `C6v`; `tools/ah_characters.py`) and the
polyhedral tables of
[Transcribe Altmann & Herzig irrep matrices for the polyhedral groups](https://github.com/hydeik/pointgroup-tables/issues/17)
(`T`, `O`, `I`; `tools/ah_matrices.py`). Their file layouts differ and are
described in the two sections below.

## Axial tables (tickets #19 and #22)

Generator matrices from the tables T n.7 for `D3`, `C4v`, `D3h`, `D2d`, `D6h`
and `C6v` (ticket #19) and for `D5`, `D8`, `D9`, `D10` (T 25.7, T 28.7,
T 29.7, T 30.7, the tables that `D5h`, `D8h`, `D10h`, `C5v`, `C8v`, `C10v`
share; ticket #22): every two-dimensional irrep, ordinary and projective, on
the book's principal rotation, its first secondary element and (where
present) `sigma_h` or `i`. Read from the page images; they verify the
spin-harmonic rule of the `altmann-herzig` basis convention
(`docs/spec/verification.md`, `REF-7`; `tools/ah_axial_rule.py check`
compares them with the rule and with the book's own kets in `../bases/`).

```
uv run python tools/ah_characters.py generate   # writes these files from the MATRICES dict of the tool
uv run python tools/ah_characters.py check      # unitarity, D(g)^order = +-1, trace = character table entry
```

The trace check reads the transcribed character table of the group; the
`D_n` tables have none (`D5`, `D8`, `D9`, `D10` are not among the 39
character files), so their traces are compared with the computed characters
of the canonical irrep of the same symbol.

### File layout

```toml
group = "C4v"
table = "T 52.7, p. 490"
status = "done"
conductor = 8
book_epsilon = "exp(2 pi i / 8)"      # the table's epsilon (e in book_values)
elements = ["C4+", "sigma_v1"]        # book labels of the generators, in order

[irreps."E1/2"]                       # the book's irrep label
book_values = [                       # one matrix per element, the book's tokens (overbar = minus)
  [["e*", "0"], ["0", "e"]],
  [["0", "-i"], ["-i", "0"]],
]
matrices = [                          # the same in the spec's value grammar for `conductor`
  [["zeta^7", "0"], ["0", "zeta"]],
  [["0", "-i"], ["-i", "0"]],
]
```

A table with several roots of unity names them all instead of
`book_epsilon`, with the tokens used in `book_values`
(`d` = delta, `e` = epsilon, `h` = eta, `t` = theta of the table's footnote):

```toml
book_symbols = { d = "exp(2 pi i / 9)", e = "exp(4 pi i / 9)", h = "exp(6 pi i / 9)", t = "exp(8 pi i / 9)" }
```

## Polyhedral groups (ticket #17)

One TOML file per rotation group whose tables the book prints: `T.toml`
(T 70.7, p. 592), `O.toml` (T 69.7, pp. 582-584), `I.toml` (T 74.7a/b,
pp. 654-655). The other polyhedral groups reuse them (book section 16-7,
p. 78): `Td` reads T 69.7 with its own operation names (T 73.7, `S4x-` for
`C4x+`, `sigma_d1` for `C2'a`, ...); `Th`, `Oh`, `Ih` are `H x Ci` and use
`D_g(i h) = D(h)`, `D_u(i h) = -D(h)` (T 72.7, T 71.7, T 75.7). Spinor
irreps are tabulated on the untilded operations, `D(~g) = -D(g)`.

```
uv run python tools/ah_matrices.py check       # verify the transcription held in the tool
uv run python tools/ah_matrices.py generate    # (re)write T.toml, O.toml, I.toml
```

### How the files were produced

The tool holds the printed entries in the book's own symbols (transcribed
from page images of the PDF, since the OCR layer scrambles matrix layout)
and evaluates them exactly in `Q(zeta_M)`. For `O` and `T` the book prints
one matrix per operation; for `I` it prints the matrices of three
generators (`C5,1+`, `C3,1+`, `C2a`) and a word in those generators for
every operation. The irrational constants of T 74.7b are printed as
`r exp(i arctan x)`; they were identified as exact surds in
`Q(sqrt2, sqrt3, sqrt5)` (PSLQ at 80 digits, residual below `1e-80`) and the
tool re-evaluates the printed formulas in floating point against them.

Before a file is written, the tool checks exactly, in the field:

- every matrix is unitary;
- the matrices multiply as `D(a) D(b) = sign * D(c)`, the sign being that of
  the quaternion product of the book's own parameters (T 71.1 / T 75.1, the
  untilded lifts of ticket #13), i.e. `-1` exactly where the book's factor
  table has a tilde, for every pair of operations (`O`, `T`) or every
  operation times every generator (`I`, matrices along the words);
- the trace of every matrix equals the character of exactly one canonical
  irrep of the same dimension, which fixes the `[canonical]` map of book
  symbols to canonical Mulliken symbols;
- every entry survives a round trip through the value grammar.

### File layout

```toml
group = "O"
table = "T 69.7, pp. 582-584"
status = "done"
conductor = 24                 # every entry lies in Q(zeta_24); "zeta" below is zeta_24
notes = """..."""

[book]
irreps = ["E", "T1", "T2", "E1/2", "E5/2", "F3/2"]   # the book's symbols, its column order
elements = ["E", "C2x", ...]   # book labels of the transcribed matrices (labels/<Group>.toml maps them)
generators = ["C4z+", "C2'd"]  # a generating set among them

[book.constants]               # the footnote of the table, as printed
delta = "exp(2 pi i/24)"

[book.exact]                   # the same constants in the value grammar; "ds" is the conjugate of "d"
d = "zeta"

[[known_discrepancy]]          # I only: a printed matrix corrected before use (docs/spec/verification.md)
element = "C5,1+"
irrep = "F3/2"
fix = "transpose"
printed = "sqrt2/4 | eta -i*sm mus -ls; ..."

[canonical]                    # book irrep symbol -> canonical Mulliken symbol (by the traces)
"F3/2" = "G3/2"

[minimal_conductor]            # smallest Q(zeta_d) that holds every printed entry of the irrep
"F3/2" = 24

[words]                        # I only: T 74.7a, each operation as a product of generators, left to right
"C5,2+" = ["C5,1+", "C5,1+", "C2a", "C3,1+"]

[irreps."E1/2"]
matrices = [                   # one per entry of [book].elements, rows then columns, value grammar
  [ ["1", "0"], ["0", "1"] ],
]
```

The value grammar is that of `docs/spec/basis-conventions.md`: `i`,
`sqrt2`, `sqrt3`, `sqrt5`, `phi` and `zeta` (= `zeta_M` for the file's
`conductor`), integers, `+ - * / ^` and parentheses. The tool writes each
entry in the shortest of: a rational, `r*zeta^k`, `r*sqrtN*zeta^k`, or
`re + i*im` with `re`, `im` rational combinations of `1`, `sqrt2`, `sqrt3`,
`sqrt5` and their products.

### Facts recorded for the basis convention

- **Untilded operations.** The book's spinor matrices are those of its
  untilded operations (`D(C4z+) = diag(epsilon*, epsilon)` in `E1/2` of
  `O`, trace `sqrt2`; `D(C5,1+)` of `E1/2` of `I` has trace `(1 + sqrt5)/2`),
  so a generator image can be lifted by the canonical section wherever the
  book's untilded operation is the canonical one (ticket #13 lists the
  exceptions per group as `swapped`).
- **Words.** The words of T 74.7a are words in the rotation group: 45 of the
  60 multiply to the tilded operation in the double group (`E = C2a C2a` is
  `~E`), so products along them give `-D(g)` for those operations in the
  spinor irreps. The generator matrices themselves are untilded.
- **Errata of T 74.7b.** With `delta = exp(2 pi i/24)` as printed, the
  `C3,1+` matrices of `F3/2` and `I5/2` cube to `+1` on their delta blocks
  (a spinor `C3` must cube to `-1`; T 70.7 prints the same blocks with the
  opposite sign), and the `C5,1+` matrices of `F3/2` and `I5/2` are printed
  with rows and columns interchanged. Reading `delta` as `-delta` in those
  `C3,1+` matrices and transposing those `C5,1+` matrices is the only
  correction (among single signs, conjugations and transpositions) under
  which every product along the words agrees with the factor system and
  every trace with the character. `I.toml` carries the corrected matrices
  and records the printed forms under `[[known_discrepancy]]`.
- **Minimal conductors** (smallest `Q(zeta_d)` per irrep): `T`: `T` 4,
  `E1/2` 4, `1F3/2` and `2F3/2` 12; `O`: `E` 3, `T1` and `T2` 4, `E1/2` and
  `E5/2` 8, `F3/2` 24; `I`: `T1`, `T2`, `F`, `E1/2`, `E7/2` 20, `F3/2` 60,
  `H` and `I5/2` 120. A convention keeping the book's entries therefore
  declares 24 for `2T` and `2O` and 120 for `2I`.
- **Symbols.** By the traces, the book's `F3/2` is the canonical `G3/2`
  (`O`, `I`), `F` of `I` is `G`, `1F3/2` and `2F3/2` of `T` are `1E3/2` and
  `2E3/2` in that order.
- **No implicit bases.** Every irrep of dimension >= 2 has printed matrices
  (or printed generator matrices); the book's own account of the bases is
  T *n*.6 (symmetrized bases / spin harmonics), which the transcription did
  not need.

## Redistribution

Committed to the repository as test data by the decision of
[Verification strategy](https://github.com/hydeik/pointgroup-tables/issues/12)
(ADR-0006). See `tests/reference/README.md` for the licensing note.
