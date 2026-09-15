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

## Amendment (ticket #22)

The transcription of every symmetrized basis the book prints
(`tests/reference/altmann-herzig/bases/`, 28 tables) and of four more
matrix tables (`D5`, `D8`, `D9`, `D10`) showed the rule as first stated,
from six sample tables, to be wrong in three places, all corrected in
`docs/spec/basis-conventions.md` and `tools/ah_axial_rule.py`:

- the chain divides `N` by its **smallest prime factor** at each step
  (`3, 9` for `n = 9`, `3, 9, 18` for `D9h`), not by two down to the odd
  part; the odd part alone gave the wrong first partner for `E5/2`,
  `E13/2` and four ordinary irreps of the `n = 9` groups;
- the phase of the secondary element is the residue at the **deciding
  level**, not at the bottom level (`E3` of `D9`, `E3'`, `E3''`, `E9/2` of
  `D9h`, `E3` of `D6d`);
- `E_{N/4}` of `D6d` and `D10d` takes the negative partner first, and the
  ordinary irreps of the odd `Cnv` groups carry the `D_n` matrices of
  `C2'(k)` although their planes lie at azimuth `π/2` (the book maps
  `sigma_v(k)` to `C2'(k)` and picks the sign of the second ket to match).

The decision stands: the amended rule reproduces every printed ket and
matrix for `n ≤ 10` exactly, and the transcription is test data that
checks it. The cost noted above materialised once (a systematic
misreading, caught by the transcription task the option anticipated); the
extrapolation beyond the book keeps two unconfirmed points, named in the
spec.
