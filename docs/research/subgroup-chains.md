# Subgroup catalogue: which group chains of every family are canonical

Research note for wayfinder ticket
[#41](https://github.com/hydeik/pointgroup-tables/issues/41), part of
[Map: version 2, chain-adapted bases and Clebsch–Gordan coefficients](https://github.com/hydeik/pointgroup-tables/issues/38).
Vocabulary as in `GLOSSARY.md` ("Group chains and coupling"): a **group chain**
`G ⊃ H ⊃ …` is **canonical** when every restriction is multiplicity-free, and
the **chain labels** of a component are the Mulliken symbols of the groups
below `G`.

Every statement below is tagged:

- **[computed]**: produced by `docs/research/scripts/subgroup_chains.py` with
  the version-1 library (characters of `G` and of `H`), exhaustively for the
  axial families at `n = 1..12` and the seven polyhedral groups. The full
  machine-generated tables (subgroup classes, branching rules, chains) are in
  `docs/research/subgroup-chains-catalogue.md` (600 KB, regenerate with
  `uv run python3 docs/research/scripts/subgroup_chains.py --out docs/research/subgroup-chains-catalogue.md`).
- **[proven]**: follows for every `n` from the character formulas; the proof
  is given or sketched.
- **[guess]**: an expectation not covered by either.

## 1. Method

**Subgroups.** The subgroups of the abstract group are enumerated from the
Cayley table (closure of the cyclic subgroups under joins with cyclic
subgroups), conjugacy classes of subgroups by conjugating element sets. A
subgroup is identified as a point group by reading the family table of
`docs/spec/abstract-groups.md` backwards: its proper elements form a rotation
group `S` (cyclic, dihedral, `T`, `O` or `I` by order and maximal element
order); no improper element gives type I, the inversion gives type II
(`Cnh`/`S2n`, `Dnh`/`Dnd` by the parity of `n`), otherwise type III (the
rotation parts of all elements form `R`, the proper ones its index-2
subgroup `H`). Degenerate symbols are normalised through the alias table
(`D1 = C2`, `D1h = C2v`, `C1v = C1h`).

**Embedding.** Each subgroup class is matched to the library's group of that
symbol by an explicit isomorphism of double groups `H_lib† → H†` (`H†` the
preimage of `H` in `G†`) that preserves parity, the geometric invariants of
every element read off its canonical label (kind, rotation order,
`|exponent|`) and, for every element that is not a half turn, which preimage
is the lifted one (half turns have a conventional section). A small
generating set of `H_lib†` (highest element order first) is sent to the
admissible images of **smallest canonical index**; the representative of a
conjugacy class is the conjugate whose generator images are
lexicographically smallest. This is deterministic and is one candidate for
the embedding rule of ticket #42 (section 7).

**Branching.** For every irrep `Γ` of `G` and every irrep `γ` of `H_lib`,
`m = (1/|H†|) Σ_{h ∈ H†} χ_Γ(φ(h)) conj(χ_γ(h))`, exactly in the cyclotomic
field of conductor `lcm(M_G, M_H)`; the sum `Σ m·dim γ = dim Γ` is checked.
A step `G ⊃ H` is **free** for a kind of irrep when every `m ≤ 1` over that
kind.

**Chains.** A chain descends through proper subgroups (not only maximal
ones), never through the trivial group `C1` (a step to `C1` is never free for
projective irreps and free for ordinary irreps only from an abelian group),
and stops at the first group whose rotation group is cyclic (`Cn`, `Cnh`,
`S2n`, including `C1h`, `S2`): its double group is abelian, so every ordinary
and projective irrep is one-dimensional and no further step can add a label.
Chains are enumerated through steps that are free for at least one kind of
irrep, and reduced to one per orbit under simultaneous conjugation in `G`. A
chain is canonical for ordinary / projective / both irreps when all its steps
are. The column "abelian from" of the catalogue names the first group of the
chain whose ordinary irreps are all one-dimensional: from there the chain
labels of ordinary irreps are complete; for projective irreps the labels are
complete only at the terminal group (`D2`, `C2v`, `D2h` are abelian but have
the two-dimensional projective irrep `E1/2`).

Sanity checks passed by the computation: the subgroup class counts of `O`
(11 including `O` and `C1`), `Oh` (33), `I` (9), `Ih` (22) agree with the
literature; the branching rules reproduce the textbook ones (`Oh ⊃ D4h`:
`T1g → A2g + Eg`, `Oh ⊃ Td`: `A1u → A2`, `T1u → T2`; `O ⊃ D3`: `T1 → A2 + E`;
`I ⊃ T`: `G → A + T`, `H → 1E + 2E + T`; `I ⊃ D5`: `H → A1 + E1 + E2`).

## 2. Two structural reductions [proven]

Everything reduces to the rotation groups `Cn`, `Dn`, `T`, `O`, `I`.

**Twisted groups are their rotation group.** A type III group
`G = {(r, χ(r)) : r ∈ R}` is mapped isomorphically onto `R` by
`(r, p) ↦ r`, and `(q, p) ↦ q` maps `G†` onto `R†` (the parity is a function
of the rotation part). Subgroup lattice, conjugacy classes, characters and
therefore every branching multiplicity of `G` are those of `R`; only the
symbols change. The dictionary:

| `G` | `R` | subgroup of `R` | subgroup of `G` |
| --- | --- | --- | --- |
| `Cnv` | `Dn` | `Cd`; `Dd` (`d ∣ n`); `D1` (a `C2'`) | `Cd`; `Cdv`; `C1h` (a `sigma_v`/`sigma_d`) |
| `Cnh`, `n` odd | `C2n` | `C2d` (`d ∣ n`); `Cd` | `Cdh` (`C1h` for `d = 1`); `Cd` |
| `S2n`, `n` even | `C2n` | `Cd` (`d ∣ n`); `C2d` (`d ∣ n`, `n/d` odd) | `Cd`; `S2d` (`S4` in `S12`, `S8` in `S24`) |
| `Dnh`, `n` odd | `D2n` | `C2n`; `C2d`, `Cd` (`d ∣ n`); `D2d` (`d ∣ n`, one class); `Dd` (`d ∣ n`, two classes) | `Cnh`; `Cdh` (`C1h`, the `sigma_h`, for `d = 1`), `Cd`; `Ddh` (`C2v` for `d = 1`); `Dd` and `Cdv` (`C2'` and `C1h` (`sigma_v`) for `d = 1`) |
| `Dnd`, `n` even | `D2n` | `C2n`; `Cd` (`d ∣ n`); `C2d` (`d ∣ n`, `n/d` odd); `D2d` (`d ∣ n`, `n/d` odd, one class); `Dd` (`d ∣ n`, two classes) | `S2n`; `Cd`; `S2d`; `Ddd`; `Dd` and `Cdv` (`C2'` and `C1h` (`sigma_d`) for `d = 1`) |
| `Td` | `O` | `T`, `D4`, `D3`, `D2`, `D2'`, `C4`, `C3`, `C2`, `C2'` | `T`, `D2d`, `C3v`, `D2`, `C2v`, `S4`, `C3`, `C2`, `C1h` |

**Centro groups are a direct product.** A type II group is
`G = R × {±1}`, `G† = R† × {±1}`, and its irreps are `Γ_g = Γ ⊗ 1`,
`Γ_u = Γ ⊗ sgn`. Its subgroups are the pairs `(S, S0)`: `S ≤ R` proper
(`H = S`), centro (`H = S × {±1}`), or twisted with an index-2 subgroup
`S0 < S` (`H ≅ S`). In each case `Γ_g ↓ H` and `Γ_u ↓ H` are `Γ ↓ S` up to a
one-dimensional character, so **the step `G ⊃ H` is free exactly when the
step `R ⊃ S` is**, for both kinds of irrep, and `H` has as many conjugates as
`(S, S0)` has in `R`. Examples [computed]: `Oh ⊃ D2h` (the normal `D2h`,
`S = D2` normal in `O`) fails for both kinds like `O ⊃ D2`; `Oh ⊃ D2h'`
(`S = D2' = ⟨C2, C2'⟩`) and `Oh ⊃ S6` (`S = C3`) are free for ordinary and
fail for projective irreps like `O ⊃ D2'` and `O ⊃ C3`; `Oh ⊃ Td` (`S = O`)
and `Oh ⊃ C4v` (`S = D4`) are free for both.

Consequences: the axial families have two shapes, the **cyclic shape** (`Cn`;
`Cnh` and `S2n` of either parity, as `C2n` or `Cn × Ci`) and the **dihedral
shape** `D_N` (`Dn` and `Cnv` with `N = n`; `Dnh` for odd `n` and `Dnd` for
even `n` with `N = 2n`; `Dnh` for even `n` and `Dnd` for odd `n` as
`Dn × Ci`). The polyhedral groups are `T`, `O`, `I` and `Td ≅ O`,
`Th = T × Ci`, `Oh = O × Ci`, `Ih = I × Ci`.

## 3. The cyclic shape: `Cn`, `Cnh`, `S2n` [proven, computed for n ≤ 12]

Every subgroup is normal (one conjugate) and is again of the cyclic shape:
`Cd` for `d ∣ n` in `Cn`; in `Cnh` (`n` even) `Cd`, `Cdh` (`d` even) and
`S2d` (`d` odd) for `d ∣ n`; in `Cnh` (`n` odd) `Cd` and `Cdh` for `d ∣ n`
(`C1h` for `d = 1`); in `S2n` with `n` odd (`= Cn × Ci`) `Cd` and `S2d` for
`d ∣ n`; in `S2n` with `n` even (`≅ C2n`) `Cd` for `d ∣ n` and `S2d` for the
`d ∣ n` with `n/d` odd (`S12 ⊃ S4`, `S24 ⊃ S8`, but no `S4` in `S8`).

Every step is free for both kinds of irrep, because every irrep (ordinary or
projective) of every group of the shape is one-dimensional and a
one-dimensional representation restricts to one irrep. The catalogue confirms
`max multiplicity = 1` for every step of `C1..C12`, `C1h..C12h`, `S2..S24`.

Chains of this shape are therefore all canonical, and they are redundant:
the labels of `Cn` already distinguish every component (the double group is
abelian). The only thing a chain decides here is the **partner orientation**:
the restriction `1E_k ↓ Cd` is `1E_{k mod d}` only if the principal element of
`Cd` is identified with `C_n^{n/d}` and not with its inverse (section 7).

## 4. The dihedral shape `D_N`

Let `N` be the axis order of the shape (`N = n` for `Dn`, `Cnv`; `N = 2n`
for `Dnh` with odd `n` and `Dnd` with even `n`; `N = n` with the product by
`Ci` for `Dnh` with even `n` and `Dnd` with odd `n`). Irreps of `D_N`:
ordinary `A1`, `A2`, `E_k` (`1 ≤ k < N/2`), plus `B1`, `B2` for even `N`;
projective `E_{j/2}` for odd `j < N` (two-dimensional), plus the
one-dimensional pair `1E_{N/2}`, `2E_{N/2}` for odd `N`.

### 4.1 Subgroup classes [proven; computed for N ≤ 12 and the derived families]

| subgroup of `D_N` | condition | classes | conjugates per class | distinguishing canonical elements |
| --- | --- | --- | --- | --- |
| `C_d` | `d ∣ N` | 1 (normal) | 1 | `C_N^{N/d}` (the principal element `C_N^1` for `d = N`) |
| `D_d` | `d ∣ N`, `N/d` odd | 1 | `N/d` | the class contains both `C2'(1)` and `C2''(1)` |
| `D_d` | `d ∣ N`, `N/d` even | 2 | `N/(2d)` each | one class contains `C2'(1)`, the other `C2''(1)` (for `d = 1`: the `C2'` and the `C2''` subgroups) |

`D_1 = C2` is the subgroup generated by one secondary element (`C2'`, or a
reflection in the derived families); `D_2` is the abelian `D2`. The counts are
read off the dihedral group: `D_d = ⟨C_N^{N/d}, C_N^k C2'⟩` depends on
`k mod N/d`, and conjugation by `C_N` shifts `k` by 2, so the classes are the
residues of `k mod gcd(2, N/d)`. The catalogue lists the classes of every
derived family with their symbols (e.g. `D6h`: 31 proper classes, `D12h`:
61, `D6d`: 14), all consistent with this table through section 2.

### 4.2 Single steps [proven; computed for N ≤ 12]

| step | ordinary irreps | projective irreps | both |
| --- | --- | --- | --- |
| `D_N ⊃ D_d`, any `d ∣ N`, `d < N` (either class; includes `D_N ⊃ C2'`) | free | free | **free** |
| `D_N ⊃ C_N` | free | free | **free** |
| `D_N ⊃ C_d`, `d ∣ N`, `d < N`, `d` odd, `N = 2d` | free | multiplicity 2 | no |
| `D_N ⊃ C_d`, `d ∣ N`, `d < N`, `d` odd, `N ≥ 3d` | multiplicity 2 | multiplicity 2 | no |
| `D_N ⊃ C_d`, `d ∣ N`, `d < N`, `d` even (includes the central `C2 = ⟨C_N^{N/2}⟩`) | multiplicity 2 | free | no |

Proofs (characters on the generator `c = C_N^{N/d}` of `C_d`):

- `E_k ↓ C_d = (k) + (−k) mod d`; the two characters coincide iff `d ∣ 2k`.
  For `d = N` no `k < N/2` qualifies: **free**. For `d < N`: if `d` is even,
  `k = d/2 < N/2` qualifies; if `d` is odd, `k = d` qualifies iff `d < N/2`,
  i.e. `N/d ≥ 3`. The one-dimensional irreps restrict to one irrep each.
- `E_{j/2} ↓ C_d†`: the two components have eigenvalues `exp(±iπ j/d)` on the
  lifted `c`, equal iff `d ∣ j`. `j` is odd, so this needs `d` odd; `d = N`
  is excluded by `j < N`; for odd `d < N` the irrep `E_{d/2}` exists
  (`d ≤ N/3 < N − 2` when `N` is odd) and fails. Hence **free iff `d = N` or
  `d` even**. The one-dimensional pair of odd `N` restricts to the single
  irrep `B_{N/2}` of `C_N†` (character `−1` on the lifted `C_N`): still
  multiplicity-free, two irreps of `G` sharing one label of `H`.
- `D_N ⊃ D_d`: restricted to `C_d`, `E_k` is either an irreducible
  `E_{k mod d}` of `D_d` or, when `d ∣ 2k`, the sum of two one-dimensional
  irreps that differ in the sign on `C2'` (the character of `E_k` on every
  secondary element is 0): `A1 + A2` or `B1 + B2`, distinct. Likewise
  `E_{j/2}` restricts to one `E_{j'/2}` or to `1E_{d/2} + 2E_{d/2}`.
  **Free for both, always.**

So for `Dn ⊃ Cn` (the ticket's example): free for ordinary irreps for every
`n`, and free for projective irreps for every `n`; the parity of `n` only
decides whether all projective irreps are two-dimensional (even `n`) or two
of them are one-dimensional and share the `Cn` label `B_{n/2}` (odd `n`).

### 4.3 The same facts in Schoenflies symbols [computed, n ≤ 12; proven through section 2]

| family | free for both | free for ordinary only | free for projective only | fails for both |
| --- | --- | --- | --- | --- |
| `Dn` | `Cn`; every `Dd` (`d ∣ n`, `d < n`, both classes); `C2'`, `C2''` | `Cd` with `d` odd, `n = 2d` (`D6 ⊃ C3`, `D10 ⊃ C5`) | `Cd` with `d` even `< n` (`D4 ⊃ C2`, `D8 ⊃ C4`, `D12 ⊃ C6`) | `Cd` with `d` odd, `n ≥ 3d` (`D9 ⊃ C3`, `D12 ⊃ C3`) |
| `Cnv` | `Cn`; every `Cdv`; `C1h` (`sigma_v`, `sigma_d`) | `Cd`, `d` odd, `n = 2d` | `Cd`, `d` even `< n` | `Cd`, `d` odd, `n ≥ 3d` |
| `Dnh`, `n` odd (`≅ D2n`) | `Cnh`; `Dd`, `Cdv`, `Ddh` (`d ∣ n`, `d < n`); `C2v`; `C2'`; `C1h` (`sigma_v`) | `Cn` (`= C_d` of `D2n` with `2n = 2d`: `D3h ⊃ C3`, `D5h ⊃ C5`) | `C1h` (`sigma_h`, the central `C2` of `D2n`); `Cdh` for `d ∣ n`, `d < n` (`C2d` of `D2n`, even: `D9h ⊃ C3h`) | `Cd` for `d ∣ n`, `d < n` (`D9h ⊃ C3`) |
| `Dnd`, `n` even (`≅ D2n`) | `S2n`; `Dd`, `Cdv` (`d ∣ n`); `Ddd` (`d ∣ n`, `n/d` odd); `C2'`; `C1h` (`sigma_d`) | none (`C_d` of `D2n` with `d` odd and `2n = 2d` would need `n` odd) | `Cd` for even `d ∣ n`, `d < 2n` (`Cn` itself: `D4d ⊃ C4`, `D6d ⊃ C6`, `D4d ⊃ C2`); `S2d` for `d ∣ n`, `n/d` odd, `d < n` (`D6d ⊃ S4`) | `Cd` for odd `d ∣ n` (`D6d ⊃ C3`, `D10d ⊃ C5`) |
| `Dnh`, `n` even (`= Dn × Ci`); the status of `H` is that of `Dn ⊃ S`, `S` the rotation group of `H` | `S = Dn`: `Dn`, `Cnv`, `D_{n/2}h` or `D_{n/2}d`; `S = Cn`: `Cn`, `Cnh`, `S_n`/`C_{n/2}h`; `S = Dd` (`d ∣ n`, `d < n`): `Dd`, `Cdv`, `Ddh`, `Ddd`, …; `S = C2'`: `C2'`, `C2h` (`C2'`, `sigma`), `C1h` (`sigma_v`/`sigma_d`); also `C2v` (`C2'`, `sigma_h`) with `S = D2'` | `S = Cd`, `d` odd, `n = 2d`: `Cd`, `S2d` (`D6h ⊃ C3`, `D6h ⊃ S6`, `D10h ⊃ C5`, `D10h ⊃ S10`) | `S = Cd`, `d` even, `d < n`: `Cd`, `Cdh`, `S_d`/`C_{d/2}h` (`D4h ⊃ C2`, `C2h`, `C1h (sigma_h)`; `D8h ⊃ C4`, `C4h`, `S4`; `D12h ⊃ C6`, `C6h`, `C3h`) | `S = C1`: `S2` (`n ≥ 4`; `D2h ⊃ S2` is free for ordinary irreps, `D2` being abelian); `S = Cd`, `d` odd, `n ≥ 3d`: `Cd`, `S2d` (`D12h ⊃ C3`, `D12h ⊃ S6`) |
| `Dnd`, `n` odd (`= Dn × Ci`) | `S = Dn`: `Dn`, `Cnv`; `S = Cn`: `Cn`, `S2n`; `S = Dd` (`d ∣ n`, `d < n`): `Dd`, `Cdv`, `Ddd`; `S = C2'`: `C2'`, `C2h`, `C1h` (`sigma_d`) | none (`n` odd has no `d` with `n = 2d`) | none (`n` odd has no even `d`) | `S = C1`: `S2` (`D3d ⊃ S2`); `S = Cd`, `d ∣ n`, `d < n`: `Cd`, `S2d` (`D9d ⊃ C3`, `D9d ⊃ S6`) |

### 4.4 Chains of the dihedral shape [proven from 4.2; computed lists for n ≤ 12 in the catalogue]

A chain of `D_N` is a descent through dihedral subgroups
`D_N ⊃ D_{d1} ⊃ D_{d2} ⊃ …` (each `d_{i+1} ∣ d_i`, either class at even
quotients) that ends with one cyclic step `D_d ⊃ C_d`, or at `D_1 = C2`
(`C2'`-type, the real `Dn ⊃ C2x` basis of Chen et al.), or at `D_2 = D2`
followed by `D2 ⊃ C2` (any of its three `C2`).

- **Canonical for both kinds** iff every cyclic step is `D_d ⊃ C_d` (same
  order). Examples [computed]: `D6 ⊃ C6`, `D6 ⊃ D3 ⊃ C3`, `D6 ⊃ D3 ⊃ C2`,
  `D6 ⊃ D2 ⊃ C2`, `D6 ⊃ C2'`; `D8 ⊃ D4 ⊃ D2 ⊃ C2`; `C6v ⊃ C3v ⊃ C3`;
  `D3h ⊃ C3h`, `D3h ⊃ D3 ⊃ C3`, `D3h ⊃ C3v ⊃ C3`, `D3h ⊃ C2v ⊃ C2`;
  `D4d ⊃ S8`, `D4d ⊃ D4 ⊃ C4`; `D4h ⊃ C4h`, `D4h ⊃ D2d ⊃ S4`,
  `D4h ⊃ D2h ⊃ C2h`; `D6h ⊃ C6h`, `D6h ⊃ D3d ⊃ S6`, `D6h ⊃ D3h ⊃ C3h`.
- **Canonical for ordinary irreps only**: a cyclic step `D_{2d} ⊃ C_d` with
  `d` odd (`D6 ⊃ C3`, `C6v ⊃ C3`, `D3h ⊃ C3`, `D5h ⊃ C5`, `D6h ⊃ C3`,
  `D6h ⊃ S6`, `D6d ⊃ D6 ⊃ C3`), and the step `D2h ⊃ S2` (`E1/2g ↓ Ci† = 2`).
- **Canonical for projective irreps only**: a cyclic step `D_m ⊃ C_d` with
  `d` even `< m` (`D4 ⊃ C2`, `D8 ⊃ C4`, `D8 ⊃ D4 ⊃ C2`, `C4v ⊃ C2`,
  `D3h ⊃ C1h(sigma_h)`, `D4h ⊃ C2h`, `D4d ⊃ C4`, `D6d ⊃ S4`).
- Where labels are complete: ordinary irreps are fully labelled once the
  chain reaches an abelian group (`D2`, `C2v`, `D2h`, `C2h`, or the cyclic
  terminal); projective irreps need the cyclic terminal (`D2†`, `C2v†`,
  `D2h†` are non-abelian, `E1/2` is two-dimensional). Chains ending in a
  secondary `C2` (`Dn ⊃ C2x`) or a reflection (`Cnv ⊃ C1h`) are complete for
  both kinds (`C2†` and `C1h†` are cyclic of order 4).
- Symbol strings are ambiguous for the two classes of `D_d` (`N/d` even), for
  the two classes of `D_1` (`C2'` vs `C2''`, `sigma_v` vs `sigma_d`), and for
  `D_m ⊃ C2` (central vs secondary `C2`: `D4 ⊃ C2` names three chain classes
  with different canonicity). The two `D_d` classes have the same canonicity
  and isomorphic branching with `B1 ↔ B2` exchanged (`D4 ⊃ D2`:
  `B1 → A, B2 → B1` for the `C2'` class, `B1 → B1, B2 → A` for the `C2''`
  class) [computed].

## 5. The polyhedral groups [computed]

Notation: `H (c)` is a class of `c` conjugates; the representative is the
conjugate generated by the listed canonical elements of `G`; "o/p" is the
maximal multiplicity for ordinary/projective irreps (`1/1` = free for both).

### 5.1 `T` (order 12)

| class | conjugates | generators in `T` | maximal in | o/p |
| --- | --- | --- | --- | --- |
| `D2` | 1 | `C2(1)[1], C2(2)[2]` | `T` | 1/1 |
| `C3` | 4 | `C3^1(1)[4]` | `T` | 1/1 |
| `C2` | 3 | `C2(1)[1]` | `D2` | 2/1 (`T → A + 2B`) |

Maximal subgroups: `D2`, `C3`. Branching: `T ⊃ C3`: `T → A + 1E + 2E`,
`E1/2 → 1E1/2 + 2E1/2`, `1E3/2 → B3/2 + 1E1/2`, `2E3/2 → B3/2 + 2E1/2`;
`T ⊃ D2`: `T → B1 + B2 + B3`, every projective irrep `→ E1/2`.
Chains: **both** `T ⊃ C3`, `T ⊃ D2 ⊃ C2`; **projective only** `T ⊃ C2`.
`T ⊃ D2 ⊃ C2` has one class: the three `C2` of `D2` are conjugate in `T`.

### 5.2 `Td` (order 24, `≅ O`)

| class | conjugates | generators in `Td` | maximal in | o/p |
| --- | --- | --- | --- | --- |
| `T` | 1 | `C3^1(1)[4], C3^-1(2)[6]` | `Td` | 1/1 |
| `D2d` | 3 | `S4^-1(1)[18], C2(2)[2]` | `Td` | 1/1 |
| `C3v` | 4 | `C3^1(1)[4], sigma_d(2)[13]` | `Td` | 1/1 |
| `D2` | 1 | `C2(1)[1], C2(2)[2]` | `T`, `D2d` | 2/2 |
| `C2v` | 3 | `C2(1)[1], sigma_d(5)[16]` | `D2d` | 1/2 |
| `S4` | 3 | `S4^-1(1)[18]` | `D2d` | 1/1 |
| `C3` | 4 | `C3^1(1)[4]` | `T`, `C3v` | 1/2 (`G3/2 → 2B3/2 + …`) |
| `C2` | 3 | `C2(1)[1]` | `D2`, `C2v`, `S4` | 2/2 |
| `C1h` | 6 | `sigma_d(1)[12]` | `C3v`, `C2v` | 2/2 |

Chains: **both** `Td ⊃ S4`; `Td ⊃ C3v ⊃ C3`, `Td ⊃ C3v ⊃ C1h`;
`Td ⊃ D2d ⊃ S4`, `Td ⊃ D2d ⊃ C1h`, `Td ⊃ D2d ⊃ C2` (the `C2` perpendicular
to the `S4` axis), `Td ⊃ D2d ⊃ C2v ⊃ C2|C1h`, `Td ⊃ D2d ⊃ D2 ⊃ C2` (2
classes); `Td ⊃ T ⊃ C3`, `Td ⊃ T ⊃ D2 ⊃ C2`. **Ordinary only** `Td ⊃ C3`,
`Td ⊃ C2v ⊃ C2|C1h`. **Projective only** `Td ⊃ D2d ⊃ C2` (the `C2` on the
`S4` axis), `Td ⊃ T ⊃ C2`. Ambiguous strings: `Td ⊃ D2d ⊃ C2` (2 classes,
different canonicity), `Td ⊃ D2d ⊃ D2 ⊃ C2` (2 classes, both canonical).

### 5.3 `Th` (order 24, `= T × Ci`)

| class | conjugates | generators in `Th` | maximal in | o/p |
| --- | --- | --- | --- | --- |
| `T` | 1 | `C3^1(1)[4], C3^1(2)[5]` | `Th` | 1/1 |
| `D2h` | 1 | `C2(1)[1], C2(2)[2], sigma_h(1)[13]` | `Th` | 1/1 |
| `S6` | 4 | `C3^1(1)[4], S6^-1(1)[16]` | `Th` | 1/1 |
| `D2` | 1 | `C2(1)[1], C2(2)[2]` | `T`, `D2h` | 1/1 |
| `C2h` | 3 | `C2(1)[1], sigma_h(1)[13]` | `D2h` | 2/1 |
| `C2v` | 3 | `C2(1)[1], sigma_h(2)[14]` | `D2h` | 1/1 |
| `C3` | 4 | `C3^1(1)[4]` | `T`, `S6` | 1/1 |
| `C2` | 3 | `C2(1)[1]` | `D2`, `C2h`, `C2v` | 2/1 |
| `S2` | 1 | `i[12]` | `S6`, `C2h` | 3/2 |
| `C1h` | 3 | `sigma_h(1)[13]` | `C2h`, `C2v` | 2/1 |

Chains: **both** `Th ⊃ C3`, `Th ⊃ S6`; `Th ⊃ C2v ⊃ C2|C1h`,
`Th ⊃ D2 ⊃ C2`, `Th ⊃ D2h ⊃ C2|C2h|C1h`, `Th ⊃ T ⊃ C3`,
`Th ⊃ D2h ⊃ C2v ⊃ C2|C1h`, `Th ⊃ D2h ⊃ D2 ⊃ C2`, `Th ⊃ T ⊃ D2 ⊃ C2`.
**Ordinary only** `Th ⊃ D2h ⊃ S2`. **Projective only** `Th ⊃ C2`,
`Th ⊃ C2h`, `Th ⊃ C1h`, `Th ⊃ T ⊃ C2`. Ambiguous: `Th ⊃ C2v ⊃ C1h` and
`Th ⊃ D2h ⊃ C2v ⊃ C1h` (2 classes each, both canonical: the two planes of
one `C2v` are not conjugate in `Th`).

### 5.4 `O` (order 24)

| class | conjugates | generators in `O` | maximal in | o/p |
| --- | --- | --- | --- | --- |
| `T` | 1 | `C3^1(1)[10], C3^-1(2)[12]` | `O` | 1/1 |
| `D4` | 3 | `C4^1(1)[18], C2(2)[2]` | `O` | 1/1 |
| `D3` | 4 | `C3^1(1)[10], C2'(2)[5]` | `O` | 1/1 |
| `D2` (normal) | 1 | `C2(1)[1], C2(2)[2]` | `T`, `D4` | 2/2 |
| `D2'` | 3 | `C2(1)[1], C2'(5)[8]` | `D4` | 1/2 |
| `C4` | 3 | `C4^1(1)[18]` | `D4` | 1/1 |
| `C3` | 4 | `C3^1(1)[10]` | `T`, `D3` | 1/2 |
| `C2` | 3 | `C2(1)[1]` | `D2`, `D2'`, `C4` | 2/2 |
| `C2'` | 6 | `C2'(1)[4]` | `D3`, `D2'` | 2/2 |

Branching of the free steps: `O ⊃ C4`: `E → A + B`, `T1 → A + 1E + 2E`,
`T2 → B + 1E + 2E`, `E1/2 → 1E1/2 + 2E1/2`, `E5/2 → 1E3/2 + 2E3/2`,
`G3/2 → 1E1/2 + 2E1/2 + 1E3/2 + 2E3/2`; `O ⊃ T`: `E → 1E + 2E`,
`E5/2 → E1/2`, `G3/2 → 1E3/2 + 2E3/2`; `O ⊃ D3`: `T1 → A2 + E`,
`T2 → A1 + E`, `G3/2 → E1/2 + 1E3/2 + 2E3/2`; `O ⊃ C3` fails only through
`G3/2 → 2B3/2 + 1E1/2 + 2E1/2`.
Chains: **both** `O ⊃ C4`; `O ⊃ D4 ⊃ C4`, `O ⊃ D4 ⊃ C2` (`C2(2)` or
`C2'(5)`, i.e. a `C2` perpendicular to the `C4` axis), `O ⊃ D3 ⊃ C3`,
`O ⊃ D3 ⊃ C2`, `O ⊃ T ⊃ C3`; `O ⊃ D4 ⊃ D2 ⊃ C2` (4 classes),
`O ⊃ T ⊃ D2 ⊃ C2`. **Ordinary only** `O ⊃ C3`, `O ⊃ D2' ⊃ C2` (2 classes).
**Projective only** `O ⊃ D4 ⊃ C2` (the `C2` on the `C4` axis),
`O ⊃ T ⊃ C2`. Ambiguous strings: `O ⊃ D4 ⊃ C2` (3 classes: one
projective-only, two both), `O ⊃ D2 ⊃ C2`, `O ⊃ D4 ⊃ D2 ⊃ C2`. The
`O ⊃ D4 ⊃ C4` chain adds nothing to `O ⊃ C4` (every `D4` irrep restricts to
distinct `C4` irreps, so the `D4` label is a function of the `C4` label);
Ping & Chen 1999 call them equivalent.

### 5.5 `Oh` (order 48, `= O × Ci`)

32 proper classes (33 with `Oh`), maximal: `O`, `Th`, `Td`, `D4h` (3),
`D3d` (4). Free for both: `O`, `Th`, `Td`, `D4h` (3), `T`, `D3d` (4), `D4` (3),
`D2d` (both classes: on the `C2` axes with `S4^-1(1), C2(2)`, and with
`S4^-1(1), C2'(5)`), `C4h` (3), `C4v` (3), `D3` (4), `C3v` (4), `C4` (3),
`S4` (3). Free for ordinary only: `D2h'` (3, `C2(1), C2'(5), sigma_h(1)`),
`S6` (4), `D2'` (3), `C2v'` (3, `C2(1), sigma_d(5)`), `C2v''` (6,
`C2'(1), sigma_h(3)`), `C3` (4). Fail for both: the normal `D2h` and `D2`,
`C2h` (3), `C2v` (3, `C2(1), sigma_h(2)`), `C2h'` (6, `C2'(1), sigma_d(1)`),
`C2` (3), `C2'` (6), `S2`, `C1h` (3, `sigma_h`), `C1h` (6, `sigma_d`). Every
status equals that of the rotation subgroup in `O` (section 2).

Chains: 212 chain classes, 160 canonical for both, 32 for ordinary only, 20
for projective only. Canonical for both with at most two steps:
`Oh ⊃ C4`, `Oh ⊃ C4h`, `Oh ⊃ S4`; `Oh ⊃ C3v ⊃ C3|C1h`, `Oh ⊃ C4v ⊃ C4|C1h`,
`Oh ⊃ D2d ⊃ S4|C2|C1h`, `Oh ⊃ D3 ⊃ C3|C2`, `Oh ⊃ D3d ⊃ S6|C3|C2h|C2|C1h`,
`Oh ⊃ D4 ⊃ C4|C2`, `Oh ⊃ D4h ⊃ C4h|C4|S4|C2h|C2|C1h`, `Oh ⊃ O ⊃ C4`,
`Oh ⊃ T ⊃ C3`, `Oh ⊃ Td ⊃ S4`, `Oh ⊃ Th ⊃ C3|S6`; longer ones include
`Oh ⊃ O ⊃ T ⊃ C3`, `Oh ⊃ Td ⊃ T ⊃ C3`, `Oh ⊃ Th ⊃ T ⊃ C3`,
`Oh ⊃ T ⊃ D2 ⊃ C2`, `Oh ⊃ O ⊃ D4 ⊃ C4`, `Oh ⊃ D4h ⊃ D2h ⊃ C2h`,
`Oh ⊃ D4h ⊃ D2d ⊃ S4`. 49 symbol strings name more than one chain class
(e.g. `Oh ⊃ D4h ⊃ C2`: 3 classes, one projective-only; `Oh ⊃ D4h ⊃ C2v ⊃
C1h`: 6 classes), see the catalogue.

### 5.6 `I` (order 60)

| class | conjugates | generators in `I` | maximal in | o/p |
| --- | --- | --- | --- | --- |
| `T` | 5 | `C3^1(1)[16], C3^1(5)[22]` | `I` | 1/1 |
| `D5` | 6 | `C5^1(1)[36], C2(7)[7]` | `I` | 1/1 |
| `D3` | 10 | `C3^1(1)[16], C2(10)[10]` | `I` | 2/2 (`H → A1 + 2E`) |
| `C5` | 6 | `C5^1(1)[36]` | `D5` | 1/2 (`I5/2 → 2B5/2 + …`) |
| `D2` | 5 | `C2(1)[1], C2(7)[7]` | `T` | 2/3 |
| `C3` | 10 | `C3^1(1)[16]` | `T`, `D3` | 2/2 |
| `C2` | 15 | `C2(1)[1]` | `D5`, `D3`, `D2` | 3/3 |

Branching: `I ⊃ D5`: `T1 → A2 + E1`, `T2 → A2 + E2`, `G → E1 + E2`,
`H → A1 + E1 + E2`, `E7/2 → E3/2`, `G3/2 → E1/2 + E3/2`,
`I5/2 → E1/2 + E3/2 + 1E5/2 + 2E5/2`; `I ⊃ T`: `G → A + T`,
`H → 1E + 2E + T`, `I5/2 → E1/2 + 1E3/2 + 2E3/2`; `I ⊃ C5` is free for
ordinary irreps (`H → A + 1E1 + 2E1 + 1E2 + 2E2`) and fails only through
`I5/2 → 2B5/2 + 1E1/2 + 2E1/2 + 1E3/2 + 2E3/2` (six dimensions into the five
projective irreps of `C5`).
Chains: **both** `I ⊃ D5 ⊃ C5`, `I ⊃ D5 ⊃ C2`, `I ⊃ T ⊃ C3`,
`I ⊃ T ⊃ D2 ⊃ C2`. **Ordinary only** `I ⊃ C5`. **Projective only**
`I ⊃ T ⊃ C2`. No chain through `D3` is canonical.

### 5.7 `Ih` (order 120, `= I × Ci`)

20 proper classes (22 with `Ih` and `C1`), maximal: `I`, `Th` (5), `D5d` (6),
`D3d` (10). Free for both: `I`, `Th`, `D5d`, `T` (5), `D5` (6), `C5v` (6).
Ordinary only: `S10` (6), `C5` (6). Fail for both: `D3d`, `D2h` (5), `D3`,
`S6`, `C3v` (10 each), `C3` (10), `D2` (5), `C2h`, `C2v`, `C2`, `C1h` (15
each), `S2`. Representatives contain `C5^1(1)[36]`, `C3^1(1)[16]` or
`C2(1)[1]`.

Chains: **both** `Ih ⊃ C5v ⊃ C5|C1h`, `Ih ⊃ D5 ⊃ C5|C2`,
`Ih ⊃ D5d ⊃ S10|C5|C2h|C2|C1h`, `Ih ⊃ T ⊃ C3`, `Ih ⊃ Th ⊃ C3|S6`,
`Ih ⊃ D5d ⊃ C5v ⊃ C5|C1h`, `Ih ⊃ D5d ⊃ D5 ⊃ C5|C2`, `Ih ⊃ I ⊃ D5 ⊃ C5|C2`,
`Ih ⊃ I ⊃ T ⊃ C3`, `Ih ⊃ T ⊃ D2 ⊃ C2`, `Ih ⊃ Th ⊃ C2v ⊃ C2|C1h`,
`Ih ⊃ Th ⊃ D2 ⊃ C2`, `Ih ⊃ Th ⊃ D2h ⊃ C2h|C2|C1h`, `Ih ⊃ Th ⊃ T ⊃ C3`,
`Ih ⊃ I ⊃ T ⊃ D2 ⊃ C2`, `Ih ⊃ Th ⊃ D2h ⊃ C2v ⊃ C2|C1h`,
`Ih ⊃ Th ⊃ D2h ⊃ D2 ⊃ C2`, `Ih ⊃ Th ⊃ T ⊃ D2 ⊃ C2`. **Ordinary only**
`Ih ⊃ C5`, `Ih ⊃ S10`, `Ih ⊃ I ⊃ C5`, `Ih ⊃ Th ⊃ D2h ⊃ S2`. **Projective
only** `Ih ⊃ T ⊃ C2`, `Ih ⊃ Th ⊃ C2|C2h|C1h`, `Ih ⊃ I ⊃ T ⊃ C2`,
`Ih ⊃ Th ⊃ T ⊃ C2`. Ambiguous: `Ih ⊃ Th ⊃ C2v ⊃ C1h` and
`Ih ⊃ Th ⊃ D2h ⊃ C2v ⊃ C1h` (2 classes each, both canonical).

### 5.8 Maximal-subgroup lattices (classes)

```
T   : T > {D2, C3};  D2 > C2;  C3, C2 > C1
Td  : Td > {T, D2d(3), C3v(4)};  T > {D2, C3};  D2d > {D2, C2v, S4};  C3v > {C3, C1h};
      D2 > C2;  C2v > {C2, C1h};  S4 > C2
Th  : Th > {T, D2h, S6(4)};  T > {D2, C3};  D2h > {D2, C2h(3), C2v(3)};  S6 > {C3, S2};
      C2h > {C2, S2, C1h};  C2v > {C2, C1h};  D2 > C2
O   : O > {T, D4(3), D3(4)};  T > {D2, C3};  D4 > {D2, D2'(3), C4};  D3 > {C3, C2'};
      D2 > C2;  D2' > {C2, C2'};  C4 > C2
I   : I > {T(5), D5(6), D3(10)};  T > {D2, C3};  D5 > {C5, C2};  D3 > {C3, C2};  D2 > C2
Oh  : Oh > {O, Th, Td, D4h(3), D3d(4)};  then the lattice of O with each class
      of (S, S0) (section 2); the catalogue column "maximal in" lists it in full
Ih  : Ih > {I, Th(5), D5d(6), D3d(10)};  likewise from I
```

## 6. The chains of the Chen–Fan–Ping papers [computed status]

Papers under `PDF/` (titles and abstracts read with `pypdf`); "status" is
this note's result for the chain they use.

| paper | chain(s) | status |
| --- | --- | --- |
| Chen & Fan, J. Math. Phys. 39, 467 (1998) | `T ⊃ C3`, s.v. and d.v. | canonical for both |
| Chen & Fan, J. Math. Phys. 39, 5502 (1998) | `T ⊃ D2 ⊃ C2` (s.v.), `T ⊃ C2` (d.v.) | `T ⊃ D2 ⊃ C2` both; `T ⊃ C2` projective only (`T → A + 2B`), as the paper's split implies |
| Chen & Fan, J. Math. Phys. 39, 5519 (1998) | CG of `T†` in `T ⊃ C3`; also `T ⊃ D2 ⊃ C2`, `T ⊃ C2` | as above |
| Fan & Chen, Int. J. Quantum Chem. 74, 7 (1999) | `O ⊃ C4` (s.v. and d.v.) | both |
| Fan, Int. J. Quantum Chem. 83, 259 (2001) | `O ⊃ T ⊃ C3` | both |
| Fan, Int. J. Quantum Chem. 91, 541 (2003) | CG in `O† ⊃ C4†` | both |
| Chen, Fan & Paldus, Int. J. Quantum Chem. 76, 585 (2000) | `G ⊃ T ⊃ D2 ⊃ C2` for `G = O, Td, Th, Oh`; `O ⊃ D4 ⊃ D2 ⊃ C2` | all both (every class of `O ⊃ D4 ⊃ D2 ⊃ C2`) |
| Fan, Chen & Draayer, Acta Cryst. A55, 871 (1999) | `I ⊃ C5`, s.v. | ordinary only (fails for `I5/2`) |
| Fan, Chen & Draayer, Acta Cryst. A55, 1049 (1999) | `I ⊃ D5 ⊃ C5`, spinor | both; the switch from `I ⊃ C5` is forced by `I5/2 ↓ C5† ∋ 2B5/2` |
| Chen & Ping, J. Mol. Spectrosc. 186, 403 (1997) | `Ih ⊃ C5` bases for `C60` (s.v.) | ordinary only |
| Chen, Fan, McAven & Butler, J. Math. Phys. 41, 8196 (2000) | all dihedral groups `Dn`, `Cnv`, `Dnd`, `Dnh`: `Dn ⊃ Cn` (complex basis) and `Dn ⊃ C2x` (real basis), `Cnv ⊃ Cn`, `Cnv ⊃ Cs`; improper groups from `Dn` | `Dn ⊃ Cn`, `Dn ⊃ C2'`, `Cnv ⊃ Cn`, `Cnv ⊃ C1h` both for every `n`; the improper analogues `Dnh ⊃ Cnh`, `Dnd ⊃ S2n` both for every `n` |
| Ping & Chen, Int. J. Quantum Chem. 75, 67 (1999) | any canonical chain of the 32 crystallographic groups (program PG); Table XVI: `D3 ⊃ C3`, `D3 ⊃ C2`, `C3v ⊃ C3`, `C3v ⊃ Cs`, `D4 ⊃ D2` (s.v.), `D4 ⊃ C2x`, `C4v ⊃ C2v` (s.v.), `C4v ⊃ Cs`, `D2d ⊃ D2` (s.v.), `D2d ⊃ C2v` (s.v.), `D2d ⊃ C2x`, `D2d ⊃ Cs`, `T ⊃ C3`, `T ⊃ C2x` (d.v.), `T ⊃ D2` (s.v.) | all these steps are free for the kind stated; the paper notes the 12 conjugate choices of `O ⊃ D4 ⊃ C2`, which is this note's ambiguity of that string |
| Chen & Ping, Comput. Phys. Commun. 120, 71 (1999) | program for any canonical `G ⊃ H ⊃ K`, s.v. and d.v., crystallographic groups | method only |
| Ping et al., Comput. Phys. Commun. 52, 355 (1989); Fan et al., Int. J. Quantum Chem. 75, 1 (1999) | space groups; Euler angles | no chain |

**Canonical chains no paper treats** (both kinds): `O ⊃ D3 ⊃ C3` and
`O ⊃ D3 ⊃ C2` (Butler's tables use `O ⊃ D3 ⊃ C3`, not the Chen group);
`O ⊃ D4 ⊃ C2` (perpendicular `C2`); `Td ⊃ S4`, `Td ⊃ C3v ⊃ C3`,
`Td ⊃ D2d ⊃ S4`; `Th ⊃ S6`, `Th ⊃ D2h ⊃ C2h`; `Oh ⊃ C4h`, `Oh ⊃ S4`,
`Oh ⊃ D3d ⊃ S6`, `Oh ⊃ D4h ⊃ C4h`; `I ⊃ T ⊃ C3`, `I ⊃ D5 ⊃ C2`,
`I ⊃ T ⊃ D2 ⊃ C2` (Altmann & Herzig tabulate `I ⊃ T ⊃ C2`-adapted
functions); `Ih ⊃ D5d ⊃ S10`, `Ih ⊃ D5d ⊃ C5`, `Ih ⊃ Th ⊃ S6`; for the
dihedral shape every `D_N ⊃ D_d ⊃ … ⊃ C_d` descent (`D6 ⊃ D3 ⊃ C3`,
`D4h ⊃ D2d ⊃ S4`, `D6h ⊃ D3d ⊃ S6`, …) and the mixed terminals
`Dnh ⊃ Cnv ⊃ Cn`, `Dnd ⊃ Dn ⊃ Cn`.

## 7. What this means for the chain grammar and the default chain

Facts and options only; the decision is ticket #42.

**Which chains are canonical for both kinds.**

- Dihedral shape (all axial families through section 2): exactly the chains
  whose only cyclic step is `D_d ⊃ C_d` of the same order (or that end at a
  secondary `C2`/reflection). In Schoenflies terms the one-step chains
  `Dn ⊃ Cn`, `Cnv ⊃ Cn`, `Dnh ⊃ Cnh`, `Dnd ⊃ S2n` are canonical for both
  for every `n`, and the natural "divide the axis order" descents
  `D_N ⊃ D_{N/p} ⊃ …` are canonical as long as the final cyclic step keeps
  the order. Cyclic shape: everything.
- Polyhedral: `T ⊃ C3`, `T ⊃ D2 ⊃ C2`; `O ⊃ C4`, `O ⊃ T ⊃ C3`,
  `O ⊃ D4 ⊃ C4`, `O ⊃ D3 ⊃ C3`, `O ⊃ D4 ⊃ D2 ⊃ C2`, `O ⊃ T ⊃ D2 ⊃ C2`;
  `Td ⊃ S4`, `Td ⊃ T ⊃ C3`, `Td ⊃ C3v ⊃ C3`, …; `Th ⊃ S6`, `Th ⊃ C3`,
  `Th ⊃ T ⊃ C3`, …; `Oh ⊃ C4h`, `Oh ⊃ C4`, `Oh ⊃ O ⊃ C4`, `Oh ⊃ O ⊃ T ⊃ C3`,
  …; `I ⊃ D5 ⊃ C5`, `I ⊃ T ⊃ C3`, `I ⊃ T ⊃ D2 ⊃ C2`; `Ih ⊃ D5d ⊃ S10`,
  `Ih ⊃ D5d ⊃ C5`, `Ih ⊃ I ⊃ D5 ⊃ C5`, `Ih ⊃ Th ⊃ S6`, …
- Not canonical for both although a paper uses them for one kind: `T ⊃ C2`
  (projective only), `I ⊃ C5` and `Ih ⊃ C5` (ordinary only), `O ⊃ C3`
  (ordinary only), `Dn ⊃ Cd` with `d < n`.

**Where the embedding rule is non-trivial.**

1. *Conjugate subgroups.* Polyhedral: `C3`, `D3`, `C3v`, `S6`, `D3d` (4 in
   `O`-type groups, 10 in `I`-type), `C4`, `D4`, `D4h`, `C4h`, `C4v`, `S4`,
   `D2d` (3), `C5`, `D5`, `D5d`, `S10`, `C5v` (6), `T`, `Th` (5 in `I`-type),
   and all `C2`-type subgroups. Axial: `D_d` and the derived `Dd`, `Cdv`,
   `Ddh`, `Ddd` with `d < N`, every `C2'`, `C2v`, `C2h`, `C1h` built on a
   secondary element. The catalogue's representative is the conjugate
   containing the smallest-index admissible generator; for the polyhedral
   groups that is the conjugate containing the polyhedral principal element
   `C4^1(1)`, `C3^1(1)`, `S4^-1(1)` or `C5^1(1)` of its order, or `C2(1)`,
   `C2'(1)`, `sigma_h(1)`, `sigma_d(1)`; for the dihedral shape the one
   containing `C2'(1)` (option: or `sigma_v(1)`/`sigma_d(1)`, which has the
   smallest axis index but not always the smallest element index, e.g.
   `D3h ⊃ C3v` picks `sigma_v(3)[9]` by index).
2. *Non-conjugate subgroups with the same symbol.* `D_d` of the dihedral shape
   at even `N/d` (the `C2'` and the `C2''` class), `C2'` vs `C2''`,
   `sigma_v` vs `sigma_d`, and the central vs secondary `C2`; in `O`/`Oh`
   the normal `D2`/`D2h` vs `D2'`/`D2h'`, `C2` vs `C2'`, `C2v` vs `C2v'`
   vs `C2v''`, `sigma_h` vs `sigma_d`. A Schoenflies string then names
   several chain classes, sometimes with different canonicity (`D4 ⊃ C2`,
   `O ⊃ D4 ⊃ C2`, `Td ⊃ D2d ⊃ C2`, `Oh ⊃ D4h ⊃ C2`). Options: (a) a fixed
   rule picks one class per string (e.g. the smallest-index generator, which
   picks the secondary `C2'(1)` for `Dn ⊃ C2`, a both-canonical choice, and
   the central `C2(1)` for `O ⊃ D4 ⊃ C2`, a projective-only choice) and the
   other classes are unreachable; (b) the grammar admits primes or axis
   indices (`D4 ⊃ C2'`, `O ⊃ D4 ⊃ C2''`, `Dnh ⊃ C1h(sigma_h)`); (c) the
   step's canonicity decides when the classes differ (prefer the canonical
   one). In the dihedral shape the two `D_d` classes give isomorphic
   branching with `B1 ↔ B2` exchanged, so option (a) is harmless there.
3. *Orientation inside a unique subgroup.* Even a normal subgroup needs its
   principal element identified with an element of `G`: `C_n^1` vs `C_n^-1`
   decides `1E` vs `2E` labels (the cyclic shape, `Cnh ⊃ Cn`, `Dn ⊃ Cn`), and
   for a dihedral `H` which `C2` is `H`'s principal element and which is
   `C2'(1)` decides `B1`/`B2`/`B3` (`D4 ⊃ D2`, `T ⊃ D2`). The smallest-index
   rule inherits `G`'s principal element and axis numbering wherever
   possible (`C_n^1`, then `C2'(1)`).

**Default chain, options the data supports.**

- One rule for all axial families: `G ⊃ (the maximal subgroup of the cyclic
  shape with the same `N`)`: `Dn ⊃ Cn`, `Cnv ⊃ Cn`, `Dnh ⊃ Cnh`,
  `Dnd ⊃ S2n`; the cyclic shape has the trivial chain (no subgroup needed,
  or `Cnh ⊃ Cn`). All canonical for both, all one step, all normal
  subgroups, all the chains of Chen et al. (2000). Alternative: the real
  `Dn ⊃ C2'` basis, also canonical for both.
- Polyhedral, one step where it exists: `T ⊃ C3`, `O ⊃ C4`, `Td ⊃ S4`,
  `Th ⊃ S6`, `Oh ⊃ C4h`; for `I` and `Ih` no one-step chain is canonical
  for both (`I ⊃ C5` fails for `I5/2`), the shortest are `I ⊃ D5 ⊃ C5` (the
  papers' spinor chain) and `Ih ⊃ D5d ⊃ S10` / `Ih ⊃ D5d ⊃ C5`.
- Polyhedral, following the papers: `T ⊃ C3`, `O ⊃ C4` (or
  `O ⊃ T ⊃ C3`), `I ⊃ D5 ⊃ C5`, and for the centro/twisted groups the chain
  of their rotation group prefixed by `G ⊃ R`-type steps (`Oh ⊃ O ⊃ C4`,
  `Th ⊃ T ⊃ C3`, `Td ⊃ T ⊃ C3`, `Ih ⊃ I ⊃ D5 ⊃ C5`), or the uniform
  `G ⊃ T ⊃ D2 ⊃ C2` of Chen, Fan & Paldus for the cubic groups.
- Complete one-dimensional labels for both kinds need a cyclic-shape
  terminal; `T ⊃ D2`, `O ⊃ D4 ⊃ D2` stop at a two-dimensional `E1/2`.

## 8. Reproducing

```
uv run python3 docs/research/scripts/subgroup_chains.py --groups O D6h      # one or more groups to stdout
uv run python3 docs/research/scripts/subgroup_chains.py --summary           # single-step table of the catalogue
uv run python3 docs/research/scripts/subgroup_chains.py --out docs/research/subgroup-chains-catalogue.md
```

The catalogue (axial families `n = 1..12`, seven polyhedral groups) takes
about 15 s. Element indices in it are the canonical indices of the
`elements` table (`uv run pointgroup-tables elements O --double-group`).
