# The canonical basis versus the chain-adapted basis through the principal element

Research note for wayfinder ticket
[#49](https://github.com/hydeik/pointgroup-tables/issues/49), part of
[Map: version 2, chain-adapted bases and Clebsch–Gordan coefficients](https://github.com/hydeik/pointgroup-tables/issues/38).
Written 2026-10-28. Vocabulary follows `GLOSSARY.md`; the chain-adapted
basis is the one defined in [Group chains](../spec/group-chains.md), sections
"Chain-adapted basis" and "Exact computation" (tickets #43 and #44,
ADR-0011 and ADR-0013); the canonical basis is that of
[Basis conventions](../spec/basis-conventions.md) (ADR-0004).

**Question.** For which groups does the version-1 canonical basis coincide
with the chain-adapted basis of the chain through the principal element,
exactly or up to phases? Where they differ, by which diagonal unitary, or
by more? And can the default chain's basis be made identical to the
canonical basis without changing version 1?

Every claim below is **[observed]** (an exact computation of
`docs/research/scripts/canonical_versus_chain.py`, raw output in
`canonical-versus-chain-survey.txt`) unless marked **[proven]** (a short
argument given here) or **[inferred]** (a consequence drawn from the
computation and the spec, not itself computed).

---

## Summary

1. **One-step chains to the principal element's cyclic-shape group: the
   vectors are the canonical vectors, with phase `1` and `N = 1`, for every
   group of the survey** (`Dn ⊃ Cn`, `Cnv ⊃ Cn`, `Dnh ⊃ Cnh`, `Dnd ⊃ S2n`
   for `n = 2 … 12`, `T ⊃ C3`, `Th ⊃ S6`, `Td ⊃ S4`, `O ⊃ C4`, `Oh ⊃ C4h`;
   ordinary and projective irreps). No designated component has `N ≠ 1`
   and no diagonal unitary other than the identity appears **[observed]**.
   The "Chains through the principal element" property of the spec holds
   for the vectors.
2. **The component order agrees with the canonical order for `Dn`, `Cnv`,
   `Dnd` at every `n` and for `Dnh` at even `n`, so `matrices` under the
   default chain equals the canonical table entry by entry there**
   (every irrep of these 39 chains) **[observed]**. It does **not** agree for
   **`Dnh` with odd `n`**: the two components of `E_k'` with `k` odd and of
   `E_k''` with `k` even are swapped (`E'` of `D3h`; `E1'`, `E2''` of `D5h`;
   `E1'`, `E3'`, `E2''` of `D7h`; …). The spec's sentence "for the axial
   families the component order agrees too" is **wrong for odd `Dnh`**: the
   canonical basis orders components by the exponent of the principal
   element `c = S_n^-1` (order `2n`), the chain labels `1E_k'`/`2E_k'` of
   `Cnh` by the character at the proper principal element `c₀ = C_n`, and
   the two orders differ exactly when the leading exponent at `c` is `n − k`
   rather than `k` (§3.2, **[observed `n ≤ 11`, proven for every odd `n`]**).
   Consequently the default-chain document of an odd `Dnh` is **not**
   byte-identical to the version-1 document, contrary to ADR-0013's
   consequence ("default-chain documents byte-identical to version 1"):
   the matrices of the swapped irreps are conjugated by the swap.
3. **Polyhedral one-step defaults: same vectors, different order.** The
   permutation (chain component → canonical component) is `(1, 0, 2)` on
   `T` of `T` (`Tg`, `Tu` of `Th`), `(1, 0)` on `E`, `(1, 0, 2)` on `T1` and
   `(1, 2, 0, 3)` on `G3/2` of `O` (and `Td`, `Oh` `g`/`u`); every other irrep
   is in canonical order already, so its matrices are the canonical ones
   entry by entry **[observed]**.
4. **The longer chains `O ⊃ D4 ⊃ C4`, `Oh ⊃ D4h ⊃ C4h`, `Td ⊃ D2d ⊃ S4`
   span the canonical label spaces and differ by a diagonal unitary on
   the projective irreps only**: `λ = (1, −i)` on every `E1/2` block of
   `D4` and `(1, i)` on every `E3/2` block (`E1/2`, `E5/2`, `G3/2` of `O`);
   the ordinary irreps have phase `1`. Every `N = 1`. The phase is forced by
   the within-block rule: the canonical `D4` writes `D(C2'(1))` of `E1/2` as
   `[[0, i], [i, 0]]`, the canonical `O` writes `D(C2(2))` as
   `[[0, −1], [1, 0]]` **[observed]**.
5. **`I ⊃ D5 ⊃ C5` and `Ih ⊃ D5d ⊃ S10` differ by phases on every irrep
   but `A`, and by a genuine `2 × 2` mixing on `I5/2`.** The second
   component of every `E1`, `E2`, `E1/2`, `E3/2` block of `D5` carries a
   phase (`−ζ_20^4` on `E1 ⊂ T1`, `−ζ_20^8` on `E2 ⊂ T2`, `−ζ_20^3` on
   `E1/2 ⊂ E1/2`, …, §3.4); the designated components of these blocks have
   `N = 1`. On `I5/2` the two `−1` eigenvectors of `C5` (canonical
   components `0` and `5`) are split by `D5` into `1E5/2` and `2E5/2` with
   `N = 1/2`, `μ = (1 + i)/2`, and
   `w = (1 − i)/2 · (e_0 ∓ ζ_20 e_5)` **[observed]**. This is exactly the
   `(M, N) = (20, 1/2)` row of ADR-0013.
6. **Conductors.** Every norm factor the survey needed lies in the chain
   field `Q(ζ_L)`, `L = M_G`: `√1 = 1` everywhere, and for `I5/2`
   `N = 1/2` with `√2 ∉ Q(ζ_20)` (proven by the exhaustive sign search) and
   `(1 + i)/2 ∈ Q(ζ_20)`. The ADR-0013 rows `O>C4` / `O>D4>C4`: `8`,
   `I>D5>C5`: `20`, axial: `M_G` are confirmed for the chains surveyed
   **[observed]**.
7. **Answer to question 3 (§4.3).** Version 1 cannot change, so: the
   *vectors* of every one-step default chain are already the canonical
   ones; the *component order* can be made canonical for every one-step
   default (odd `Dnh` and the five polyhedral groups included) only by
   changing the component-order rule of #43 for cyclic-shape steps (order
   the step's labels by the exponent at the embedded principal element,
   descending, instead of the subgroup's canonical irrep order), which is
   a decision for the map; the default chains of `I` and `Ih` can never
   coincide with the canonical basis, because the within-block rule ties
   them to the canonical `D5` / `D5d` matrices, which the canonical `I`
   matrices do not restrict to. A spec statement can promise exactly the
   facts of items 1–5.

---

## 1. Setting and method

### 1.1 The two bases

The **canonical basis** of an irrep `Γ` of dimension `d ≥ 2`
([Basis conventions](../spec/basis-conventions.md)) diagonalises the
principal element `c` with the exponents descending (`ζ_N` for an ordinary
irrep, `ζ_2N` for a projective one, `N` the order of `c`), the remaining
phases fixed by the axial formula or by the frozen polyhedral generator
data of `src/pointgroup_tables/data/conventions/basis/canonical.toml`.

The **chain-adapted basis** for a resolved chain `G ⊃ H_1 ⊃ … ⊃ H_m`
([Group chains](../spec/group-chains.md#chain-adapted-basis)) is defined by
the label spaces `V_s = range(P^{σ_1}_1 ⋯ P^{σ_m}_m)` (products of the
exact isotypic projectors of the steps), the lexicographic component order
(each `σ_i` by its index in the canonical irrep order of `H_i,lib`), and the
norm-factor phase rule applied recursively: within a block of `σ_1 = σ` the
vectors make the restriction to `H_1` equal the chain-adapted matrices of
`σ` for the pulled-back chain, and the first vector of a block is
`w_t = P_t e_j / μ_F((P_t)_jj)` with `j` the designated canonical component.

The chains of this note are the chains **through the principal element**:
the default chains of every group (the maximal cyclic-shape subgroup
containing `c`, reached by the shortest canonical route) and, for the
polyhedral groups, the longer routes `O ⊃ D4 ⊃ C4`, `Oh ⊃ D4h ⊃ C4h`,
`Td ⊃ D2d ⊃ S4` that the eigenfunction-method papers use
([#39](https://github.com/hydeik/pointgroup-tables/issues/39)).

### 1.2 The script

`docs/research/scripts/canonical_versus_chain.py` (standard library and
`pointgroup_tables` only; `uv run python3 docs/research/scripts/canonical_versus_chain.py [-v] [G>H1>H2 …]`)
implements the spec's definition literally, in the chain field
`Q(ζ_L)`, `L = lcm(M_G, M_{H_i})`:

- **Embedding by rule** (`resolve_step`): every step is resolved as
  [Embedding](../spec/group-chains.md#embedding-which-subgroup-a-symbol-picks)
  prescribes: candidate images of a generating set of `H_lib` among the
  elements of the previous subgroup with the same label invariants (kind,
  rotation order, modulus of the exponent, parity), both lifts for half
  turns and the lifted preimage for every other element, extended along
  words to an injective homomorphism of double groups, checked for the
  invariants on every element, filtered by canonicity for both kinds
  (restriction multiplicities from the double-group characters), the
  lexicographically smallest image taken. The resolved embeddings are
  printed in the survey file (`C4^1 → C4^1(1)`, `C2'(1) → C2(2)` for `D4`
  in `O`; `C5^1 → C5^1(1)`, `C2'(1) → C2(7)` for `D5` in `I`;
  `S4^-1 → S4^-1(1)` for `S4` in `Td`; …).
- **Label spaces**: products of isotypic projectors summed over both
  preimages, exact rank (`0` or `1`) and the rank sum (`= d`) asserted.
- **The basis**: the designated vector `P_t e_j / μ(N)` with `μ(N)` the
  exact real square root when `N` is a square in the chain field (sign
  search with exact verification, as in the cyclotomic-closure script) and
  `(1 + i)/2` for `N = 1/2` otherwise; any other `N` would stop the script
  (none occurred). The other vectors of a block by the matrix-element
  projector built from the sub-chain realisation, which the script
  computes recursively by the same function in the library model of `H_1`
  (so `O ⊃ D4 ⊃ C4` uses the chain-adapted `D4 ⊃ C4`, itself verified equal
  to the canonical `D4`). Unit norm, membership in `V_s` and the
  **within-block equality** (step-1 restriction `=` sub-chain realisation,
  check `CB-2`) are asserted on every vector and block.
- **Comparison**: for every `w_s`, the support in canonical coordinates
  (one component: `π(s)` and the phase `λ_s`; several: the mixing vector),
  and the relation of the chain-adapted matrices to the canonical ones on
  all of `G†`: `equal` (entry by entry), `permuted` (equal after
  conjugation by `π`), `phased` (the diagonal unitary matters), `mixed`.

Running time: about 25 s for the whole survey (44 axial chains at
`n = 2 … 12`, ten polyhedral chains; `Ih ⊃ D5d ⊃ S10` on 240 double-group
elements in `Q(ζ_20)` included).

`n = 1` is not run: `D1`, `C1v`, `D1h`, `D1d` resolve to the aliases `C2`,
`C1h`, `C2v`, `C2h` ([Abstract groups](../spec/abstract-groups.md#aliases)); the cyclic-shape ones
have the empty default chain, `C2v ⊃ C2` is in the `n = 2` row.

---

## 2. Results: the axial families

### 2.1 `Dn ⊃ Cn`, `Cnv ⊃ Cn`, `Dnd ⊃ S2n` (`n = 2 … 12`), `Dnh ⊃ Cnh` (even `n`)

Every irrep, ordinary and projective: the chain labels are `(1E_k, 2E_k)`
(or `(1E_{k/2}, 2E_{k/2})`) in canonical component order, every `w_s = e_s`
exactly (phase `1`, `N = 1`), and the chain-adapted matrices equal the
canonical matrices **entry by entry** on every element of `G†`
**[observed]**. The survey file shows `matrices equal <count>, permuted 0,
phased 0, mixed 0; N in {1}` for each of these 39 chains. This is the
`CB-4` claim of [Verification](../spec/verification.md) and holds as stated
for these families.

### 2.2 `Dnh ⊃ Cnh`, odd `n`

The vectors are canonical (`w_s = e_{π(s)}`, phase `1`, `N = 1`) but the
order differs on half of the two-dimensional ordinary irreps
**[observed, `n = 3, 5, 7, 9, 11`]**:

| Group | Irreps with the two components swapped | Irreps in canonical order |
| ----- | -------------------------------------- | ------------------------- |
| `D3h` | `E'` | `E''`, `E1/2`, `E3/2`, `E5/2` |
| `D5h` | `E1'`, `E2''` | `E2'`, `E1''`, every `E_{k/2}` |
| `D7h` | `E1'`, `E3'`, `E2''` | `E2'`, `E1''`, `E3''`, every `E_{k/2}` |
| `D9h` | `E1'`, `E3'`, `E2''`, `E4''` | `E2'`, `E4'`, `E1''`, `E3''`, every `E_{k/2}` |
| `D11h` | `E1'`, `E3'`, `E5'`, `E2''`, `E4''` | the rest |

Rule: **`E_k'` with `k` odd and `E_k''` with `k` even are swapped**; the
projective irreps never are.

Why **[proven]**: in `Dnh` with odd `n` the principal element is
`c = S_n^-1` of order `N = 2n` and the proper principal element is
`c₀ = c² = C_n` ([Irrep labels](../spec/irrep-labels.md#named-elements)).
The canonical basis of a two-dimensional `E` orders its components by the
exponent `e ∈ (−n, n]` of `D(c) = diag(ζ_2n^{e}, ζ_2n^{−e})`, `e > 0`
first. The pair labels of `Cnh` (rule 8 of the irrep-label rules) are
assigned by the character at `c₀`: `1E_k` has `χ(c₀) = ζ_n^{+k}`. Since
`ζ_2n^{2e} = ζ_n^{e}`, the first canonical component carries the `Cnh`
label with exponent `e mod n` reduced into `(−n/2, n/2]`: positive, hence
`1E`, iff `e < n/2`. The exponents observed are `e = k` or `e = n − k`
(σ_h fixes the parity of `e`: `E_k'` has even `e`, `E_k''` odd `e`), so the
first component is `2E` exactly when `e = n − k`, i.e. for `E_k'` with `k`
odd and `E_k''` with `k` even. The projective irreps `E_{k/2}` have
`e = k` at `s(c)` directly (rule 9 uses `c`, not `c₀`), so their order
always agrees.

Consequence for the spec **[inferred]**: the sentence in
[Properties](../spec/group-chains.md#properties) "for the axial families
the component order agrees too, so `matrices` under the default chain
equals the canonical table entry by entry" and the `CB-4` row
("… for every axial family (`D6 --chain C6`, `D6h --chain C6h`, …)") are
false for `D3h`, `D5h`, `D7h`, …; so is the ADR-0013 consequence that the
default-chain document of every axial group is byte-identical to the
version-1 document: for odd `Dnh` the swapped irreps' matrices are
conjugated by the swap (the *values* are the same set, the *positions*
differ). The field is unaffected (`N = 1`).

---

## 3. Results: the polyhedral groups

### 3.1 One-step default chains `T ⊃ C3`, `Th ⊃ S6`, `Td ⊃ S4`, `O ⊃ C4`, `Oh ⊃ C4h`

Every `w_s` is a canonical vector with phase `1` and `N = 1`
**[observed]**. The component permutations `π` (chain component `s` ↦
canonical index; chain labels listed in chain order):

| Group | Irrep | Chain labels (chain order) | `π` | Canonical order reads |
| ----- | ----- | -------------------------- | --- | --------------------- |
| `T` | `T` | `A, 1E, 2E` | `(1, 0, 2)` | `1E, A, 2E` |
| `Th` | `Tg`, `Tu` | `Ag, 1Eg, 2Eg` (`u` alike) | `(1, 0, 2)` | `1Eg, Ag, 2Eg` |
| `O`, `Td` | `E` | `A, B` | `(1, 0)` | `B, A` |
| `O`, `Td` | `T1` | `A, 1E, 2E` | `(1, 0, 2)` | `1E, A, 2E` |
| `O`, `Td` | `G3/2` | `1E1/2, 2E1/2, 1E3/2, 2E3/2` | `(1, 2, 0, 3)` | `1E3/2, 1E1/2, 2E1/2, 2E3/2` |
| `Oh` | `Eg`, `Eu`, `T1g`, `T1u`, `G3/2g`, `G3/2u` | as `O` with `g`/`u` | as `O` | |

All other irreps (`A1`, `A2`, `T2`, `E1/2`, `E5/2` of `O`; `A`, `1E`, `2E`,
`E1/2`, `1E3/2`, `2E3/2` of `T`; their `Th`, `Td`, `Oh` counterparts) are in
canonical order and their chain-adapted matrices are the canonical
matrices entry by entry. The canonical order is the exponent of `c`
descending (`T1` of `O`: `ζ_4^{1}, 1, ζ_4^{−1}` = `1E, A, 2E`), the chain
order the canonical irrep order of `C4` (`A, B, 1E, 2E`), which is why the
two differ exactly on the irreps containing `A` or `B` together with an
`E` pair, or `B` before `A`.

### 3.2 `O ⊃ D4 ⊃ C4`, `Oh ⊃ D4h ⊃ C4h`, `Td ⊃ D2d ⊃ S4`

The label spaces are the canonical lines (every `V_s` is spanned by one
`e_k`, so every designated component has `N = 1`); the within-block rule
adds a phase on the **projective** irreps only **[observed]**:

| Irrep of `O` (`Td`; `Oh` `g`/`u`) | Chain labels | `π` | `λ` | Matrices |
| --- | --- | --- | --- | --- |
| `E` | `A1>A, B1>B` | `(1, 0)` | `(1, 1)` | permuted |
| `T1` | `A2>A, E>1E, E>2E` | `(1, 0, 2)` | `(1, 1, 1)` | permuted |
| `T2` | `B2>B, E>1E, E>2E` | id | `(1, 1, 1)` | equal |
| `E1/2` | `E1/2>1E1/2, E1/2>2E1/2` | id | `(1, −i)` | phased |
| `E5/2` | `E3/2>1E3/2, E3/2>2E3/2` | id | `(1, i)` | phased |
| `G3/2` | `E1/2>1E1/2, E1/2>2E1/2, E3/2>1E3/2, E3/2>2E3/2` | `(1, 2, 0, 3)` | `(1, −i, 1, i)` | phased |

(`A1`, `A2` trivially equal.) `Td ⊃ D2d ⊃ S4` and both parities of
`Oh ⊃ D4h ⊃ C4h` give the same table with the same phases, as the twisted
and centro inheritance predicts.

Source of the phase: the second vector of a block is
`w_s = P^σ_{s t} w_t`, which makes the restriction to `D4` equal the
canonical `D4` matrices. The canonical `D4` has
`D_{E1/2}(C2'(1)) = [[0, i], [i, 0]]` and `D_{E3/2}(C2'(1)) = [[0, i], [i, 0]]`,
whereas the canonical `O` restricted to the image `C2(2)` has
`D_{E1/2}(C2(2)) = [[0, −1], [1, 0]]` and `D_{E5/2}(C2(2)) = [[0, 1], [−1, 0]]`;
on the ordinary `E` the two agree (`[[0, 1], [1, 0]]` against
`[[0, 1], [1, 0]]` inside `T1`). The diagonal unitary `diag(1, ∓i)`
transports one into the other **[observed]**. The chain-adapted
`O ⊃ D4 ⊃ C4` matrices therefore are **not** the canonical `O` matrices on
the spinor irreps, even though the label spaces are identical; the
`O ⊃ C4` one-step chain gives the canonical ones. This is consistent with
the spec's property paragraph ("Longer chains through the principal
element … differ from the canonical vectors only by the within-block
phases of step 2").

### 3.3 `I ⊃ D5 ⊃ C5` and `Ih ⊃ D5d ⊃ S10`

Embedding resolved: `D5` through `C5^1(1)` with `C2'(1) ↦ C2(7)`, `C5`
through `C5^1(1)`. Results (`ζ = ζ_20`; `Ih` has the same table for `g` and
`u`) **[observed]**:

| Irrep | Chain labels (chain order) | `π` | `λ` | Matrices |
| --- | --- | --- | --- | --- |
| `A` | `A1>A` | id | `1` | equal |
| `T1` | `A2>A, E1>1E1, E1>2E1` | `(1, 0, 2)` | `(1, 1, −ζ^4)` | phased |
| `T2` | `A2>A, E2>1E2, E2>2E2` | `(1, 0, 2)` | `(1, 1, −ζ^8)` | phased |
| `G` | `E1>1E1, E1>2E1, E2>1E2, E2>2E2` | `(1, 2, 0, 3)` | `(1, −ζ^4, 1, −ζ^8)` | phased |
| `H` | `A1>A, E1>1E1, E1>2E1, E2>1E2, E2>2E2` | `(2, 1, 3, 0, 4)` | `(1, 1, ζ^4, 1, ζ^8)` | phased |
| `E1/2` | `E1/2>1E1/2, E1/2>2E1/2` | id | `(1, −ζ^3)` | phased |
| `E7/2` | `E3/2>1E3/2, E3/2>2E3/2` | id | `(1, ζ)` | phased |
| `G3/2` | `E1/2>…, E3/2>…` | `(1, 2, 0, 3)` | `(1, −ζ^7, 1, ζ^7)` | phased |
| `I5/2` | `E1/2>1E1/2, E1/2>2E1/2, E3/2>1E3/2, E3/2>2E3/2, 1E5/2>B5/2, 2E5/2>B5/2` | `(2, 3, 1, 4, –, –)` | `(1, ζ^3, 1, i, mixed, mixed)` | mixed |

Every designated component of an `E` block has `N = 1`. The two
one-dimensional blocks of `I5/2` have designated component `j = 0` with
`N = 1/2`, `μ = (1 + i)/2` (no `√2` in `Q(ζ_20)`), and

```
w_{1E5/2>B5/2} = ((1 − i)/2) · (e_0 − ζ_20 e_5),
w_{2E5/2>B5/2} = ((1 − i)/2) · (e_0 + ζ_20 e_5),
```

the `D5` split of the two `C5`-eigenvalue-`−1` canonical vectors (the
canonical `I5/2` has `D(C5^1(1)) = diag(−1, ζ^6, ζ^2, −ζ^8, −ζ^4, −1)`,
exponents `5, 3, 1, −1, −3, −5` of `ζ_10`), fixed by
`D_{1E5/2}(C2'(1)) = i`, `D_{2E5/2}(C2'(1)) = −i` in `D5` against
`D_{I5/2}(C2(7))_{0,5} = −ζ^4`, `_{5,0} = −ζ^6` in `I`.

The step-2 phases come from the same mechanism as in `O`: the canonical
`I` restricted to `C2(7)` is not the canonical `D5` on the `E` blocks
(`D_{T1}(C2(7))` has off-diagonal entries `ζ^6` and `−ζ^4` where the
canonical `D5` `E1` has `1` and `1`), and the within-block rule re-phases
the second component to make them equal. No one-step chain of `I` is
canonical (`I ⊃ C5` fails on `I5/2`, ticket #41), so **no chain through
the principal element of `I` or `Ih` reproduces the canonical vectors**
**[inferred from #41 and the computation]**.

### 3.4 Conductors and norm factors

Every chain of the survey ran entirely in `Q(ζ_{M_G})`: the chain fields
`L = lcm(M_G, M_{H_i})` equal `M_G` (printed as `L=` in the survey), every
`N` is `1` except the two `1/2` of `I5/2` (`I5/2g`, `I5/2u`), and the factor
`(1 + i)/2` lies in `Q(ζ_20)`. The ADR-0013 rows `O`: `C4`, `D4>C4` → `8`,
`I`: `D5>C5` → `20`, and "axial: `M_chain = M_G`" are confirmed for the
chains surveyed; the `(20, 1/2) ↦ (1 + i)/2` entry of the norm-factor
table is the one used **[observed]**. The real-square-root decision was
the floating-point reconstruction of the earlier script (with exact
verification), adequate for `N ∈ {1, 1/2}`; it proved `√(1/2) ∉ Q(ζ_20)` by
exhausting the sign patterns.

---

## 4. Answers to the ticket's items

### 4.1 Is the chain through the principal element canonical, and is the canonical basis adapted to it with the same component order?

- Every chain surveyed is canonical for both kinds (the rank check
  passed on every irrep; the embedding rule found a canonical admissible
  isomorphism for every step) **[observed]**.
- The canonical basis **is** chain-adapted (block-diagonal on the
  subgroup, each label space a canonical line) for every **one-step**
  chain through the principal element, polyhedral and axial, and for the
  longer `O`-type chains; for `I ⊃ D5 ⊃ C5` it is chain-adapted on every
  irrep except `I5/2`, where the `−1` eigenspace of `C5` is not split along
  canonical lines **[observed]**.
- The **component order** is the same for `Dn`, `Cnv`, `Dnd` (every `n`)
  and `Dnh` (even `n`); it differs for `Dnh` (odd `n`: §2.2) and for the
  polyhedral groups (§3.1, §3.2, §3.3).

### 4.2 The diagonal unitary where the bases differ only by phases, and the groups where they differ by more

- **Identity** (`λ = 1`) on every one-step chain and on the ordinary
  irreps of `O ⊃ D4 ⊃ C4`, `Oh ⊃ D4h ⊃ C4h`, `Td ⊃ D2d ⊃ S4`.
- `diag(1, −i)` on every `E1/2`-of-`D4` block and `diag(1, i)` on every
  `E3/2`-of-`D4` block of the spinor irreps of `O`, `Td`, `Oh` under the
  two-step chains.
- `I`, `Ih` under `D5 ⊃ C5` / `D5d ⊃ S10`: the phases of the table in
  §3.3 on the second component of every `E` block (eight of nine irreps
  of `I`, sixteen of eighteen of `Ih`), and **more than phases** on
  `I5/2` (`I5/2g`, `I5/2u`): the `2 × 2` block
  `((1 − i)/2) · [[1, −ζ_20], [1, ζ_20]]` on the canonical components
  `(0, 5)`.

### 4.3 Can the default chain's basis be made identical to the canonical basis without changing version 1?

Version 1 (the canonical realisations of every group, including the
subgroups' own canonical tables) is frozen; the only free things are the
rules of version 2: the component-order rule and the phase rule of #43.

**(a) Vectors.** For every one-step default chain (`Dn`, `Cnv`, `Dnh`,
`Dnd` at every `n`; `T`, `Th`, `Td`, `O`, `Oh`) the chain-adapted vectors
*are* the canonical vectors with phase `1`; nothing needs to change
**[observed]**. For `I` and `Ih` the default chain has two steps and the
within-block rule makes the `D5` blocks equal the canonical `D5` tables,
which the canonical `I` does not restrict to; the vectors differ by the
phases of §3.3 and `I5/2` mixes. Making them identical would require
either changing the canonical `I` (or `D5`) matrices, which is a version-1
change, or abandoning the within-block rule for the default chain (the
basis would then no longer contain the canonical `D5` blocks, breaking
`CB-2` and the "sub-chain block equality" property). **Not possible
under the frozen rules; a different phase rule would trade `CB-2` for it.**

**(b) Component order.** The lexicographic order by the subgroup's
canonical irrep order differs from the canonical order on the irreps
listed in §2.2 (odd `Dnh`), §3.1 (polyhedral one-step defaults). The
alternative rule that would make **every one-step default chain equal the
canonical table entry by entry** is: *at a cyclic-shape step, order the
labels by the exponent of the embedded principal element of that step
(`ζ_N` for ordinary, `ζ_2N` for projective labels, taken in the canonical
window), descending*. For `C4`, `S4`, `C4h` in `O`, `Td`, `Oh` this reads
`B, 1E, A, 2E` and `1E3/2, 1E1/2, 2E1/2, 2E3/2` (the canonical `O` order on
`E`, `T1`, `T2`, `G3/2`); for `C3`, `S6` in `T`, `Th` it reads `1E, A, 2E`; for
`Cnh` in odd `Dnh` the exponent is taken at `S_n^-1` (the group's `c`),
which is the canonical order by definition. For `I ⊃ D5 ⊃ C5` the first
step is dihedral, so the order of its blocks is untouched and the basis
still differs by the phases. The trade-off, for the map to weigh:

- *For the exponent rule*: `--chain default` of every group but `I`, `Ih`
  reproduces the version-1 `matrices` entry by entry (the ticket's
  wish, and the strongest possible statement for consumers); `CB-4`
  becomes "equal" instead of "equal up to a permutation"; the odd-`Dnh`
  discrepancy disappears without touching version 1.
- *Against it*: the order of a step's labels depends on the embedding's
  exponent assignment, not on the printed canonical irrep order of
  `H_lib` (`1E` still precedes `2E`, but `B` precedes `A`, and the
  polyhedral `A, 1E, 2E` reading of the papers is lost); the rule must be
  applied consistently inside pulled-back sub-chains (where it changes
  nothing for the chains here, since the dihedral steps are ordered by
  irrep order); it is a change to a resolved ticket (#43) and to
  `docs/spec/group-chains.md#component-order`, `CB-4` and ADR-0011.

**(c) What the spec can promise** under the present rules
**[inferred]**: "For every one-step chain to the principal element's
cyclic-shape group the chain-adapted basis vectors are the canonical basis
vectors, each with phase `1` and `N = 1`; `matrices --chain default` equals
the canonical table entry by entry for `Dn`, `Cnv`, `Dnd` (every `n`) and
`Dnh` (even `n`), equals it up to the swap of the two components of `E_k'`
(`k` odd) and `E_k''` (`k` even) for `Dnh` (odd `n`), and up to the
permutations of §3.1 for `T`, `Th`, `Td`, `O`, `Oh`. The longer chains
`O ⊃ D4 ⊃ C4`, `Oh ⊃ D4h ⊃ C4h`, `Td ⊃ D2d ⊃ S4` span the canonical label
spaces with the diagonal unitary `diag(1, ∓i)` on the `E1/2` / `E3/2`
blocks of the spinor irreps. The default chains of `I` and `Ih` span the
canonical label spaces on every irrep but `I5/2`, with the phases of §3.3,
and split the `−1` eigenspace of `C5` in `I5/2` into
`(1 − i)/2 · (e_0 ∓ ζ_20 e_5)`." Under the exponent rule of (b) the first
two clauses collapse to "entry by entry for every group but `I`, `Ih`".

### 4.4 Discrepancies with the current spec text

1. `docs/spec/group-chains.md#properties`, "Chains through the principal
   element": "For the axial families the component order agrees too, so
   `matrices` under the default chain equals the canonical table entry by
   entry" — false for `Dnh` with odd `n` (§2.2). The vector statement
   ("the chain-adapted basis vectors are the canonical basis vectors") is
   right.
2. `docs/spec/group-chains.md#component-order`: "For the axial families
   the two orders agree on every irrep (`E_k ↓ C_n = (1E_k, 2E_k)`, …)" —
   the same error; the parenthesis is right for `Cn` steps and wrong for
   the `Cnh` step of odd `Dnh`.
3. `docs/spec/verification.md`, `CB-4`: "equals the canonical table entry
   by entry for every axial family (`D6 --chain C6`, `D6h --chain C6h`, …)"
   — needs the odd-`Dnh` exception (or the rule change of §4.3 (b)).
4. ADR-0013 / `docs/spec/group-chains.md#the-conductor-of-a-chain-document`,
   "The default-chain document of an axial group therefore has the
   version-1 field, and since its matrices equal the canonical table entry
   by entry, the same exact values byte for byte" — the field claim holds,
   the byte-for-byte claim fails for odd `Dnh` (positions of the values
   differ in the swapped irreps).
5. The `(20, 1/2)` row of the norm-factor table and the "Longer chains
   through the principal element … differ … only by the within-block
   phases of step 2 and, for `I5/2`, by the `D5` split" sentence are
   confirmed; the `I5/2` split should be read as "phases on every `E`
   block **and** the split", since the phases are not confined to `I5/2`.

---

## 5. What the computation did not cover

- Axial `n ≥ 13` (the proof of §2.2 covers every odd `n` for the order
  statement; `N = 1` for every one-step chain is the cyclotomic-closure
  note's §2.6).
- Chains not through the principal element (`Dn ⊃ C2'`, `O ⊃ T ⊃ C3`,
  `T ⊃ D2 ⊃ C2`, `I ⊃ T ⊃ …`): outside the ticket; their norms are in the
  cyclotomic-closure note.
- The exact modular real-square-root procedure of the spec: the script's
  float reconstruction with exact verification is enough for `N ∈ {1, 1/2}`.
- The `altmann-herzig` basis convention (a `--chain` with a non-canonical
  `--basis-convention` is a usage error by the map's standing decision).
