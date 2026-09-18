# Altmann & Herzig character tables

Transcription of the character tables T n.4 of S. L. Altmann and P. Herzig,
*Point-Group Theory Tables*, 2nd ed., for wayfinder ticket
[Transcribe Altmann & Herzig character tables and sample axial irrep matrices](https://github.com/hydeik/pointgroup-tables/issues/19).
One TOML file per group, named by the normalised Schoenflies symbol (`Cs` is
`C1h.toml`, `Ci` is `S2.toml`); 73 groups: the 32 crystallographic groups,
`I`, `Ih`, `C5v`, `D5h`, `D5d`, `S8`, `D4d` (ticket #19) and, for
[Reference irrep symbols of the `altmann-herzig` convention beyond the transcribed character tables](https://github.com/hydeik/pointgroup-tables/issues/36),
every other axial group with `n <= 10` (`C5` `C7` `C8` `C9` `C10`, `S10` …
`S20`, `D5` `D7` `D8` `D9` `D10`, `D7h` … `D10h`, `D6d` … `D10d`, `C7v` …
`C10v`, `C5h` `C7h` … `C10h`): the book's whole range. Generated and checked by

```
uv run python tools/ah_characters.py generate [OCR_DIR]   # typeset tables from the OCR text, scanned ones from the tool's BOOK dict
uv run python tools/ah_characters.py check                # structural comparison with the computed tables (REF-1)
```

The typeset tables (41 files; `S18`, `S20`, `C9h`, `C10h` continue on a
second page) were parsed from the OCR text of the PDF; the scanned group
pages (32 files) were read from the page images and typed into
`tools/ah_characters.py`. Every file says which in its header. Verify
against the print.

## File layout

```toml
group = "D3h"
table = "T 32.4, p. 250"
status = "done"
source = "..."
transcribed_by = "..."
notes = ""
conductor = 12                      # field of `characters`; zeta = exp(2 pi i / conductor)
book_epsilon = "exp(2 pi i / 3)"    # the book's epsilon, when the table uses it (book_delta, book_eta,
                                    # book_theta, book_zeta likewise; "exp(2 pi i * 3 / 40)" with a numerator)

[[classes]]                         # the book's columns, in the book's order
label = "2S3"                       # the book's spelling, ASCII
size = 2
order = 6                           # order of the elements, from the label
parity = -1                         # +1 proper, -1 improper

[[irreps]]                          # the book's rows, in the book's order
label = "E1/2"                      # as printed (E1/2,g  1E3/2  F3/2  A3/2 ...)
projective = true
book_values = ["2", "1", "0", "0", "sqrt3", "0"]   # the book's notation: e = epsilon, d = delta, h = eta, t = theta,
                                                   # z = zeta, 2c5^2 = 2 cos(2 pi / 5)
characters = ["2", "1", "0", "0", "sqrt3", "0"]    # the same values in the spec's value grammar for `conductor`
tau = "c"                           # the book's reality column: a real, b complex, c quaternionic
```

For a projective irrep the entries are the characters of the book's untilded
elements; the tilded elements carry the opposite sign. The structural
comparison therefore matches projective rows up to one sign per class.

## Known discrepancy

`D5h.toml` carries two `[[known_discrepancy]]` entries, applied in file
order: `swap_projective = ["2S5", "2S5^2"]` (the projective rows of T 34.4
print these two entries interchanged relative to the book's own parameter
table T 39.1 and to the master table T 39.4 of `D10h`) and
`swap_labels = ["E3/2", "E7/2"]` (the two rows carry each other's symbol;
the book's bases T 34.6 and matrices T 34.7 decide), see the `note` of
each. The checker applies the documented corrections and reports the match
separately; a `swap_labels` entry does not affect the structural
comparison.
