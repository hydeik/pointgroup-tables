# Altmann & Herzig character tables

Transcription of the character tables T n.4 of S. L. Altmann and P. Herzig,
*Point-Group Theory Tables*, 2nd ed., for wayfinder ticket
[Transcribe Altmann & Herzig character tables and sample axial irrep matrices](https://github.com/hydeik/pointgroup-tables/issues/19).
One TOML file per group, named by the normalised Schoenflies symbol (`Cs` is
`C1h.toml`, `Ci` is `S2.toml`); 39 groups: the 32 crystallographic groups,
`I`, `Ih`, `C5v`, `D5h`, `D5d`, `S8`, `D4d`. Generated and checked by

```
uv run python tools/ah_characters.py generate [OCR_DIR]   # typeset tables from the OCR text, scanned ones from the tool's BOOK dict
uv run python tools/ah_characters.py check                # structural comparison with the computed tables (REF-1)
```

The typeset tables of the master groups (14 files) were parsed from the OCR
text of the PDF; the scanned group pages (25 files) were read from the page
images and typed into `tools/ah_characters.py`. Every file says which in its
header. Verify against the print.

## File layout

```toml
group = "D3h"
table = "T 32.4, p. 250"
status = "done"
source = "..."
transcribed_by = "..."
notes = ""
conductor = 12                      # field of `characters`; zeta = exp(2 pi i / conductor)
book_epsilon = "exp(2 pi i / 3)"    # the book's epsilon, when the table uses it (book_delta likewise)

[[classes]]                         # the book's columns, in the book's order
label = "2S3"                       # the book's spelling, ASCII
size = 2
order = 6                           # order of the elements, from the label
parity = -1                         # +1 proper, -1 improper

[[irreps]]                          # the book's rows, in the book's order
label = "E1/2"                      # as printed (E1/2,g  1E3/2  F3/2  A3/2 ...)
projective = true
book_values = ["2", "1", "0", "0", "sqrt3", "0"]   # the book's notation: e = epsilon, d = delta, 2c5^2 = 2 cos(2 pi / 5)
characters = ["2", "1", "0", "0", "sqrt3", "0"]    # the same values in the spec's value grammar for `conductor`
tau = "c"                           # the book's reality column: a real, b complex, c quaternionic
```

For a projective irrep the entries are the characters of the book's untilded
elements; the tilded elements carry the opposite sign. The structural
comparison therefore matches projective rows up to one sign per class.

## Known discrepancy

`D5h.toml` carries a `[[known_discrepancy]]` entry: the projective rows of
T 34.4 print the `2S5` and `2S5^2` entries interchanged relative to the
book's own parameter table T 39.1 and to the master table T 39.4 (`D10h`),
see the `note`. The checker applies the documented correction and reports
the match separately.
