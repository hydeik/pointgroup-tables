---
status: accepted
---

# The canonical basis is adapted to the principal element, with polyhedral phases frozen as shipped data

Irrep matrices are unique only up to unitary equivalence, so the tool needs
a default basis that depends on nothing but the abstract group. The
canonical basis diagonalises the principal element `c` of the family with
eigenvalue exponents descending (positive exponent first), which fixes every
basis up to one phase per vector because the eigenvalues of `c` are distinct
in every irrep. For the axial families the phases follow from one fixed
quaternion-to-SU(2) map (`a ↦ diag(ζ^k, ζ^{−k})`, `b ↦ [[0, 1], [1, 0]]`
ordinary, `[[0, i], [i, 0]]` projective). For the polyhedral families the
phases are chosen once to keep every matrix in the smallest possible
cyclotomic field and shipped as data, verified mechanically. Reference
conventions are separate data files of generator images. See
`docs/spec/basis-conventions.md`.

## Considered Options

- **Real orthogonal bases wherever the irrep is real** (Cartesian `x, y, z`
  and `d`-orbital bases, as in chemistry tables). Rejected: it conflicts
  with the irrep naming, which reads the eigenvalues of `c`, and makes
  subduction to the principal axis illegible; the complex `|j m⟩`-type
  basis is also what the primary reference uses for spinors.
- **Phases fixed by a formula** (e.g. "the first non-zero entry of each
  column of the second generator is real and positive"). Reproducible
  without data, but the normalisation introduces square roots that push
  `Ih†` from `Q(ζ_20)` to `Q(ζ_60)`–`Q(ζ_120)`, a 10–60× cost in exact
  arithmetic, for no gain in determinism once the data is frozen and
  checked.
- **Altmann & Herzig as the canonical basis.** Rejected: the canonical layer
  must be convention-free and defined for every `n`; the reference basis
  becomes one data file among others.
- **Intertwiners as the representation of reference conventions.** Storing
  `U` per irrep is compact, but books print matrices, not intertwiners;
  generator images are what a transcriber can copy, and `U` is derived for
  verification.

## Consequences

- The canonical matrices of every axial group are given by formulas for
  every `n`; the polyhedral ones are data produced and verified by a
  prototype, and any later change to that data is a breaking change of the
  output.
- One loader and one verification path serve the canonical basis and every
  reference convention.
- A reference convention that does not cover a group is an error, never a
  silent fallback.
- Reference conventions for the axial families are limited in version 1 to
  a constant `2 × 2` transform per irrep kind; an `n`- or `k`-dependent rule
  needs a vocabulary extension.
