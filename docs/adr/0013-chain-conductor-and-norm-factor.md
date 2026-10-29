---
status: accepted
---

# Chain documents declare a per-chain minimal conductor, and the canonical norm factor is a closed form or shipped data, never a search

The chain-adapted basis divides every designated vector by a canonical
norm factor `μ` with `μ · conj(μ) = N` (ADR-0011). Two things had to be
fixed for that to be exact and reproducible: which cyclotomic field a
`--chain` document declares, and which `μ` is "the" canonical one when the
real square root of `N` is not in the field. An exact survey run for ticket
#44 showed that a natural rule, "the simplest `μ = u + i·v` with `u`, `v`
short real sums of roots of unity", resolves every Clebsch–Gordan norm and
every norm of the chains of `T` and `O`, but leaves four norms of the
chains of `I` through `T` without a hit, while factors with twelve to
sixteen non-zero coefficients exist; any enumeration rich enough to reach
them is too large to be the definition. We decided: the **conductor of a
`--chain` document is the smallest multiple of the chain field
`lcm(M_G, M_{H_i})` dividing `lcm(·, 24)` in which every squared norm is a
Hermitian norm**, which is the version-1 `M_G` for every axial chain and
a per-chain-class value for the polyhedral groups (a shipped table, printed
by `list chains`); and the **canonical norm factor is the real positive
root when it exists, `(1 + i)/2` for `N = 1/2` without `√2`, and otherwise
a shipped entry `(M, N) ↦ μ`**, verified exactly, in the spirit of the
frozen polyhedral bases of ADR-0004. The exact real square root is a
modular decision procedure (sign search over residues modulo a prime
`p ≡ 1 (mod M)`, exact verification), not a floating-point one. `cg`
coefficients in a chain-adapted or reference basis are the unitary
transform of the canonical-basis ones, so no squared norm beyond the
canonical basis arises. See `docs/spec/group-chains.md`, "Exact
computation", and `docs/spec/values.md`.

## Considered Options

- **Per-chain minimal conductor; closed forms plus shipped norm factors;
  `cg` by transform** (chosen). Keeps the version-1 principle that the
  declared field is minimal (`VL-3` holds unchanged), gives every default
  chain of the axial families and of `T`, `O`, `Td`, `Th`, `Oh` the
  version-1 field (the axial default chains byte for byte, tickets #49
  and #57), halves the
  field of the `I` chains through `T` (`60`, not `120`), and makes `μ`
  depend on `(M, N)` alone with no search in the specification. The cost:
  a conductor table and a norm-factor data file to generate once (ticket
  #50) and verify, and one more place where the polyhedral groups are
  data rather than formula.
- **One conductor per group: `lcm(M_G, M_H, 24)`, i.e. `24` for `T`- and
  `O`-type groups and `120` for `I`-type groups** (the research
  recommendation). One line, no table. Rejected: `O --chain C4` would
  declare `24` for values that lie in `Q(ζ_8)`, breaking the minimality
  principle and `VL-3`, doubling or quadrupling the exact value size of
  every default-chain document of a polyhedral group, and the `I` chains
  through `T` would carry degree `32` instead of `16`.
- **Canonical norm factor by a finite enumeration** ("the simplest hit
  among `u + i·v` with short real sums"). Deterministic and
  implementation-independent, and sufficient for every Clebsch–Gordan norm
  and every chain of `T` and `O`. Rejected: insufficient for the chains of
  `I` through `T` (four norms without a hit up to four cosine terms, with
  `2.8 × 10⁷` candidates at the next stage), and every non-square case
  ends in a tie between symmetric candidates that only a lexicographic
  order breaks, giving arbitrary-looking phases. The enumeration survives
  as the documented provenance of the printed entries.
- **Canonical norm factor by the research search** (products of short sums
  solved as an `F_2` system over quadratic characters). Finds every factor
  met. Rejected as a definition: its result depends on the processing
  order and random seeds, so it would have to be frozen anyway; freezing
  the result is this decision.
- **Floating-point reconstruction of the square root with exact
  verification** (the research script). Rejected: it misses genuine
  squares with large coefficient denominators (a square of `Q(ζ_120)` with
  denominator `810000` was missed), so a "not a square" answer would not
  be a proof; the modular reconstruction is exact and no slower.
- **`cg` in a chain-adapted basis by re-running the coupling rule on that
  basis's matrices** (the wording of ticket #46). Keeps the designated
  coefficient real and positive in every basis. Rejected: the squared
  norms of the re-run are new values that no survey covered and that the
  norm-factor data would have to include for thousands of pairs, and the
  coupled copies would differ between bases although they span the same
  states; the transform needs no new factor and keeps the copies and
  exchange signs basis-independent.
