---
status: accepted
---

# Reference axial bases are a closed-form rule anchored by the label convention, not transcribed data

Altmann & Herzig write every two-dimensional irrep of an axial group in a
basis of two spherical or spin harmonics `⟨|j m⟩, ±|j −m⟩|` (T n.6) and
print the Wigner matrices of the operations in that basis (T n.7). On the
transcribed samples these bases are not a constant transform of the
canonical ones: which partner comes first and the sign of the secondary
element vary with the axis order, the irrep and the anchored plane. The
kets of T n.6 turn out to follow a closed rule in `n` and `k` (the chain
`D_n ⊃ D_{n/2} ⊃ … ⊃ D_{n₀}` decides the first partner, and the residue
at the bottom of the chain the phase of the secondary element), and the
Wigner matrices of those kets have a closed form whose only free parameter
is the azimuth of the secondary element the label convention anchors
(`0`, or `π/2` under the letter-preserving twist of ADR-0007).

We ship the axial part of the `altmann-herzig` basis convention as this
rule (`[axial] rule = "spin-harmonic"`, `anchors = "altmann-herzig"`),
defined for every `n`, verified against the book's tables for `n ≤ 10`,
and we let a basis convention depend on a label convention's anchoring for
that one parameter. See `docs/spec/basis-conventions.md`.

## Considered Options

- **Closed-form rule** (chosen). One paragraph of spec, no data file,
  every `n` covered, and the book's tables become test data that check the
  rule rather than inputs the tool needs. The costs: the rule is the book's
  construction re-derived, so a misreading would be systematic (mitigated
  by exact comparison with every transcribed table, and by the
  transcription task that reads the remaining T n.6 and T n.7 tables), and
  the basis convention has to name a label convention.
- **Explicit data** transcribed from T n.7 for `n ≤ 10` (the version 1
  vocabulary extended by an axial table keyed by family, `n`, `k`).
  Rejected: hundreds of matrices to read from page images for a result the
  rule reproduces exactly, coverage stopping at the book's `n ≤ 10`, and
  no account of *why* the entries are what they are.
- **Leave the axial families uncovered** (the state after ticket #18).
  Rejected: the axial groups are most of the tool's range, and the primary
  reference does tabulate them.
- **Duplicate the anchoring in the basis file** instead of naming the
  label convention. Rejected: the azimuth is a consequence of the
  isomorphism ADR-0007 fixed, and a second copy would drift from it; the
  tool reads the named label convention for `b` even when the user selects
  another for the labels, so `--basis-convention` alone fixes the
  matrices.
