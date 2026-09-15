---
status: accepted
---

# The exact arithmetic is the tool's own; version 1 has no runtime dependencies

ADR-0002 chose exact computation "with SymPy" over floating point. While
the specification was charted, every prototype and transcription tool did
its exact arithmetic without SymPy: a cyclotomic field `Q(ζ_M)` as rational
coefficient vectors over the power basis, with conjugation, inversion,
lifting between conductors, linear algebra (solving, nullspace, rank) and
the value grammar's printer and parser, about 360 lines in total, exercised
by every verification of the polyhedral bases and of the reference
transcriptions. SymPy would have replaced about 250 of those lines with
`QQ.cyclotomic_field` and `DomainMatrix`, but it has no complex
conjugation, imports in about a second, and weighs tens of megabytes.

We promote the prototypes' arithmetic to the implementation and ship
version 1 with **no runtime dependencies**: Python 3.13 and its standard
library only. Numeric approximations use `decimal` at 40 digits and round
once to double, which keeps the correctly rounded values ADR-0002's
consumers were promised. SymPy (or mpmath) becomes a dependency only when a
feature needs it, as a decision of that feature.

## Considered Options

- **Keep SymPy** (ADR-0002 as written): shorter number layer, `evalf` for
  the numeric form, SymPy's own pretty printer. Rejected for version 1:
  the saving is a few hundred lines of already verified code, against a
  heavy dependency that still needs a hand-written conjugation and whose
  pretty printer does not canonicalise cyclotomic expressions.
- **The tool's own arithmetic** (chosen): zero dependencies, fast start,
  the value grammar as the single spelling of exact values everywhere
  (documents, convention data, reference transcriptions, the pretty form).
- **mpmath alone** for the numeric form: rejected, `decimal` does the same
  job in forty lines.

## Consequences

- `docs/spec/values.md` defines the numeric form by `decimal`, and the
  pretty form as the value grammar rather than a "SymPy-style" string.
- mpmath is a development dependency: the test suite compares the
  `decimal` evaluation with it.
- Performance for large `n` remains the known cost of ADR-0002; the
  implementation is pure Python either way.
