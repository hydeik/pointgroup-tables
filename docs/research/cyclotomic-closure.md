# Cyclotomic closure of chain-adapted bases and Clebsch–Gordan coefficients

Research note for wayfinder ticket
[#40](https://github.com/hydeik/pointgroup-tables/issues/40), part of
[Map: version 2, chain-adapted bases and Clebsch–Gordan coefficients](https://github.com/hydeik/pointgroup-tables/issues/38).
Written 2026-10-27. Vocabulary follows `GLOSSARY.md`; the value model is that
of [Values](../spec/values.md) (one cyclotomic **field** `Q(ζ_M)` per
**document**, its **conductor** `M` fixed per group and basis convention),
decided under ADR-0002 as amended by ADR-0009 (the exact arithmetic is the
tool's own).

**Question.** Can the chain-adapted bases and the Clebsch–Gordan coefficients
of every point group and double group be written with values in one
cyclotomic field per document, and which conductor does that need?

Every claim below is marked **[proven]** (a theorem with a source, or a
short proof given here), **[observed]** (an exact computation of
`docs/research/scripts/cyclotomic_closure.py` on the groups it names) or
**[conjectured]**.

---

## Summary

1. **The only new numbers are square roots of squared norms.** In the
   projection (eigenfunction) method the unnormalised chain-adapted vector or
   Clebsch–Gordan vector is `w = P e_j` for a Hermitian projector `P` with
   entries in the document field `K = Q(ζ_M)`; its designated component is
   `w_j = P_jj = ⟨w|w⟩ = N`, a totally positive element of `K`. With the
   phase rule "designated component real and positive" every value of the
   normalised vector lies in `K(√N)`; the chain-adapted matrices lie in
   `K({√(N_i N_k)})` **[proven, §2.4]**.
2. **Such a square root lies in a cyclotomic field iff `K(√N)/Q` is abelian**
   (Kronecker–Weber), which is decidable with square roots in `K` alone
   **[proven, §2.2]**. It is *not* always the case: for `T ⊃ D2 ⊃ C2` and
   `O ⊃ D3 ⊃ C3` the squared norms `(3 ± √3)/6` and `(3 ± √3)/12` of the
   spinor chain-adapted vectors have non-abelian `K(√N)`, so **no cyclotomic
   field holds the real-positive-normalised chain-adapted vectors of these
   chains**, and for `O ⊃ D3 ⊃ C3` not even the chain-adapted matrices
   **[observed, exact proof by computation, §3.3]**.
3. **The phase freedom repairs this.** A unit-modulus rescaling puts the
   vector into `F ⊇ K` iff `N` is a Hermitian norm `μ μ̄` from `F`
   (Landherr's classification of Hermitian forms, §2.3). For every squared
   norm met in the survey (Clebsch–Gordan coefficients of `C3v`, `D4h`, `D5h`,
   `D6h`, `D7h`, `D8h`, `D9d`, `T`, `Td`, `O`, `Oh`, `I`, `Ih`, ordinary and
   projective; chain-adapted vectors of nine canonical chains of `T`, `O`, `I`) such a
   `μ` was found **exactly** in the chain field, except for `T ⊃ D2 ⊃ C2` in
   `Q(ζ_12)`, where `(3 ± √3)/6` is provably not a norm and `Q(ζ_24)` is
   needed **[observed; the non-norm statement proven, §2.3]**.
4. **Clebsch–Gordan coefficients in the canonical basis** have rational
   squared norms on every group surveyed: `1`, `1/2`, `1/3`, `2/3`, … for
   `T` and `O`, and for `I` rationals with the primes `2`, `3`, `5` and (only
   through the pairs of multiplicity two) `7`; with the real-positive phase
   rule they need `Q(ζ_24)` for `T`- and `O`-type groups, `Q(ζ_840)` for
   `I`-type groups and `Q(ζ_lcm(M, 8))` for the axial families; with a
   norm-factor phase rule the version-1 conductor `M` suffices on every
   surveyed group **[observed]**. For the axial families `N ∈ {1, 1/2}` for
   every `n` **[proven, §2.6]**.
5. **Recommendation (§5).** Keep one cyclotomic field per document; no
   quadratic extension of the value model is needed. `cg` documents:
   conductor `lcm(M, 8)` for the axial families and `lcm(M, 24)` for the
   polyhedral families with the real-positive phase rule (observed closure,
   conjectured for `n` beyond the survey only through the proof of §2.6),
   or `M` itself with a norm-factor phase rule. `--chain` documents: the
   conductor must contain the characters of every subgroup of the chain
   (`lcm(M_G, M_H, …)`, proven necessary) and the phase rule **must** be a
   norm-factor rule, because the real-positive rule leaves every cyclotomic
   field on some canonical chains; the surveyed chains then fit in
   `Q(ζ_24)` (`T`, `O`) and `Q(ζ_60)` or `Q(ζ_120)` (`I`, see §3.3).

---

## 1. Setting

- `K = Q(ζ_M)` is the document field of a group `G` in the canonical basis;
  `K^+ = K ∩ R` its maximal real subfield; `K` is a CM field (totally
  imaginary quadratic extension of the totally real `K^+`) for `M ≥ 3`.
- An element `a ∈ K^+` is **totally positive** if `σ(a) > 0` for every
  embedding `σ`. A **Hermitian norm** from `F ⊇ K` is an element
  `μ μ̄ = N_{F/F^+}(μ)`, `μ ∈ F^*`; Hermitian norms are totally positive.
- The **designated component** of a vector `w = P e_j` is `w_j`.
- Two phase rules are compared throughout: the **real-positive rule**
  (the designated component of the normalised vector is real and positive,
  i.e. the normalised vector is `w/√N`) and a **norm-factor rule** (the
  normalised vector is `w/μ` for a canonical `μ ∈ F` with `μ μ̄ = N`).
- Chains are canonical (multiplicity-free) as the map decided; the
  chain-adapted vectors are the common eigenvectors of the class operators
  of the subgroups, i.e. `P e_j` for the product `P` of the isotypic
  projectors `P^σ_H = (d_σ/|H†|) Σ_{h∈H†} χ_σ(h)^* D(h)` along the chain.

## 2. What is proven

### 2.1 Realisability over cyclotomic fields (Brauer, Reynolds; Schur indices)

- **Brauer's theorem.** Every complex representation of a finite group `G`
  of exponent `e` is realisable over `Q(ζ_e)`: R. Brauer, *On the
  representation of a group of order g in the field of the g-th roots of
  unity*, Amer. J. Math. 67 (1945) 461–471, in the exponent form of Serre,
  *Linear Representations of Finite Groups*, §12.3, Theorem 24
  (the statement is quoted in that form in
  [arXiv:2409.02589](https://arxiv.org/pdf/2409.02589) and
  [arXiv:2012.10892](https://arxiv.org/pdf/2012.10892)). Equivalently the
  Schur index of every irreducible character is `1` over `Q(ζ_e)`.
- **Projective representations.** W. F. Reynolds, *Projective
  representations of finite groups in cyclotomic fields*, Illinois J. Math.
  9 (1965) 191–198 (<https://projecteuclid.org/journals/illinois-journal-of-mathematics/volume-9/issue-2/Projective-representations-of-finite-groups-in-cyclotomic-fields/10.1215/ijm/1256067877.full>)
  carries Brauer's result to projective representations; the sharpening by
  the exponent of the representation group is J. L. Alperin and T.-N. Kuo,
  *The exponent and the projective representations of a finite group*,
  Illinois J. Math. 11 (1967). For this project the projective irreps of
  `G` are ordinary irreps of the double group `G†`, so Brauer applied to
  `G†` (exponent `12`, `24`, `60` for `2T`, `2O`, `2I`; `2N` or `4N` for the
  axial families) already gives realisability over `Q(ζ_12)`, `Q(ζ_24)`,
  `Q(ζ_60)`.
- **Schur indices of the double groups.** `SL(2,3) = 2T`, the binary
  octahedral group `2O = CSU(2,3)` and `SL(2,5) = 2I` have Schur index `2`
  for their faithful (symplectic, Frobenius–Schur indicator `−1`)
  irreps over the character field (`E1/2` of `2T`: `Q`; `E1/2`, `E5/2` of
  `2O`: `Q(√2)`; `E1/2`, `E7/2` of `2I`: `Q(√5)`; `G3/2`, `I5/2` of `2I`:
  `Q`) and index `1` for all others
  ([GroupNames, SL(2,3)](https://people.maths.bris.ac.uk/~matyd/GroupNames/1/SL(2,3).html),
  [CSU(2,3)](https://people.maths.bris.ac.uk/~matyd/GroupNames/1/CSU(2,3).html),
  [SL(2,5)](https://people.maths.bris.ac.uk/~matyd/GroupNames/97/SL(2,5).html);
  the quaternion group `Q8 = 2D2` likewise,
  [Q8](https://people.maths.bris.ac.uk/~matyd/GroupNames/1/Q8.html)).
  A Schur index `2` means the irrep is realisable over no field of odd
  degree over the character field, in particular over no real field, and
  over every quadratic extension that splits the quaternion algebra; the
  version-1 conductors `12`, `8`, `20` are such splitting fields, and
  **version 1 already ships unitary realisations over `Q(ζ_12)`, `Q(ζ_8)`,
  `Q(ζ_20)`** ([Basis conventions](../spec/basis-conventions.md)). Brauer's
  theorem is therefore not the binding constraint for this ticket: the
  irreps are fine, the question is the *normalisation* of eigenvectors.
- **What Brauer does not give.** A realisation over `Q(ζ_e)` is not
  automatically unitary with respect to the standard form; see §2.3. GAP's
  RepnDecomp, the one open-source tool that unitarises representations
  over cyclotomics, does so by an `LDL^*` decomposition of the invariant
  Hermitian form (K. Hymabaccus and D. Pasechnik, *RepnDecomp: A GAP
  package for decomposing linear representations of finite groups*, JOSS
  5(50) 1835 (2020); manual
  [chapter 3](https://docs.gap-system.org/pkg/repndecomp/doc/chap3_mj.html)),
  which produces `L D L^*` with `D` diagonal: exactly the "unnormalised
  orthogonal vectors plus their squared norms" shape of §2.4, and it stays
  in the cyclotomic field precisely because it does not take the square
  roots of `D`.

### 2.2 When is a square root cyclotomic? (Kronecker–Weber)

**Theorem (Kronecker–Weber).** A finite extension of `Q` lies in a
cyclotomic field iff it is abelian over `Q`; the smallest such field has
the conductor of the extension
([Wikipedia, Kronecker–Weber theorem](https://en.wikipedia.org/wiki/Kronecker%E2%80%93Weber_theorem),
with the attribution Kronecker 1853, Weber 1886, Hilbert 1896).

**Lemma 1 (criterion used by the script).** Let `K/Q` be abelian and
`a ∈ K^*` not a square in `K`. Then `K(√a)/Q` is Galois iff `σ(a)/a ∈ K^{*2}`
for every `σ ∈ Gal(K/Q)`; writing `σ(a)/a = b_σ^2`, it is abelian iff
`σ(b_ρ) b_σ = ρ(b_σ) b_ρ` for all `σ, ρ` (independently of the signs of the
`b`). *Proof.* Galois: an extension of `σ` to `K(√a)` must send `√a` to a
root of `X^2 − σ(a)`, which lies in `K(√a)` iff `σ(a) ∈ K^{*2} ∪ a K^{*2}`;
`σ(a) ∈ K^{*2}` would make `a` a square. Abelian: the lifts
`σ̃(√a) = b_σ √a` and the generator `τ` of `Gal(K(√a)/K)` commute with
`τ` (as `b_σ ∈ K`), and `σ̃ρ̃(√a) = σ(b_ρ) b_σ √a`, `ρ̃σ̃(√a) = ρ(b_σ) b_ρ √a`;
changing the sign of one `b` flips both sides. ∎

Consequences for the radicals the ticket names: `√3 ∈ Q(ζ_12)`,
`√7 ∈ Q(ζ_28)`, `√(2 + √2) = 2 cos(π/8) ∈ Q(ζ_16)`,
`√((5 + √5)/2) = 2 cos(π/10) ∈ Q(ζ_20)`, `√(5 + 2√5) = tan(2π/5) ∈ Q(ζ_20)`;
but `Q(√(1 + √3))` is not Galois since `(1 + √3)(1 − √3) = −2` is not a
square in `Q(√3)`, so **`√(1 + √3)` lies in no cyclotomic field**. The same
test shows `√((3 + √3)/6) = cos(θ/2)` with `cos θ = 1/√3` (half the angle
between a three-fold and a two-fold axis of the tetrahedron) is not
cyclotomic: `(3 + √3)(3 − √3) = 6` is not a square in `Q(√3)`.

A totally positive `a ∈ K^+` has a square root in `K` iff it has a *real*
one (a purely imaginary root would make `a` negative), so "square in `K`"
is decided by the exhaustive sign search of `real_sqrt` in the script,
which is a proof in both directions (§3.0).

### 2.3 Unitary realisations over a CM field: Hermitian norms

A set of pairwise orthogonal vectors `w_1, …, w_d ∈ K^d` with squared
norms `N_i` becomes an orthonormal basis after rescaling `w_i/μ_i` with
`|μ_i|^2 = N_i`. The entries of the rescaled vectors lie in a field
`F ⊇ K` iff each `μ_i` can be chosen in `F`:

**Proposition 2.** For a totally positive `N ∈ K^+` and a cyclotomic (CM)
field `F ⊇ K`: (a) `w/√N ∈ F^d` iff `√N ∈ F`; (b) there is a unit-modulus
`e^{iφ}` with `e^{iφ} w/√N ∈ F^d` iff `N ∈ N_{F/F^+}(F^*)` (`N` is a
Hermitian norm from `F`). *Proof.* (b) `e^{iφ}/√N = 1/μ` with `μ μ̄ = N`,
and `w ∈ K^d ⊆ F^d` has a non-zero component, so `1/μ ∈ F`. ∎

This is the rank-one case of the classification of Hermitian forms over a
quadratic extension `E/K_0` of number fields: two Hermitian spaces are
isometric iff they are so at every place (Landherr's local–global
principle; M. Kirschmer, *Definite quadratic and hermitian forms with
small class number*, Habilitationsschrift RWTH Aachen 2016,
[Theorem 2.4.1](http://www.math.rwth-aachen.de/~Markus.Kirschmer/papers/herm.pdf),
with the proof referred to W. Scharlau, *Quadratic and Hermitian Forms*,
Grundlehren 270, Springer 1985, Chapter 10), a Hermitian space over a
local quadratic extension is `⟨1, …, 1, d⟩` and determined by its rank and
its determinant in `K_0^*/N(E^*)` (Kirschmer, Theorem 3.2.3), and over `R`
by rank and signature (Kirschmer, Proposition 2.1.5). Hence a totally
positive definite Hermitian form over `F` is equivalent to the identity
form iff its determinant is a Hermitian norm from `F`; for the *diagonal*
rescaling we need each `N_i` is a norm individually.

**Deciding Hermitian norms.** `F/F^+` is cyclic of degree two, so by the
**Hasse norm theorem** (Hasse 1931; special case of degree two Hilbert
1897; [Wikipedia](https://en.wikipedia.org/wiki/Hasse_norm_theorem)) `N` is
a norm iff it is a local norm at every place of `F^+`. At the real places
this is total positivity; at a finite prime `P` of `F^+` unramified in `F`
it is "`ord_P(N)` even or `P` split"; ramified primes need a Hilbert
symbol. For `F = Q(ζ_M)` with `M` **not a prime power** (`M ≢ 2 mod 4`)
there are no finite ramified primes: `1 − ζ_M` is a unit (L. C. Washington,
*Introduction to Cyclotomic Fields*, 2nd ed., GTM 83, Proposition 2.8, as
cited in [arXiv:2407.02002](https://arxiv.org/html/2407.02002)), and
`ζ_M` has minimal polynomial `X^2 − (ζ + ζ^{-1})X + 1` over `F^+` with
discriminant `(ζ − ζ^{-1})^2 = −ζ^{-2}(1 − ζ^2)^2`, a unit (for `M` odd or
`4 | M` the order of `ζ^2` is again not a prime power), so the relative
discriminant is trivial. All document fields of this project except the
prime-power conductors `M ∈ {8, 16, 32, …}` (`D4h`, `O`, `Oh`, `C8`, …; `M = 2^k`)
and `M = p^k` (`C_p`, `C_p h` with odd `p`, `M = 5` for `C5`) are of this
kind.

*Worked example.* `N = (3 − √3)/6 ∈ Q(ζ_12)^+ = Q(√3)`. `3 − √3 = √3(√3 − 1)`,
`6 = (√3 − 1)^2 (√3)^2 × unit`, so `N = unit / ((√3 − 1)(√3))` has odd
valuation at the primes `P_2 = (√3 − 1)` and `P_3 = (√3)` of `Q(√3)`. In
`Q(ζ_12)` both `2` and `3` have a single prime of residue degree `2`
(`e = f = 2`), so `P_2`, `P_3` (residue degree `1`) are inert in
`Q(ζ_12)/Q(√3)`: `N` is **not** a Hermitian norm from `Q(ζ_12)`, and no
phase puts the `T ⊃ D2 ⊃ C2` spinor vector into `Q(ζ_12)^2`. In `Q(ζ_24)`
the script exhibits `μ` with `μ μ̄ = N` (§3.3), so `Q(ζ_24)` works. ∎

The script does not implement the local theory; it *exhibits* `μ` by a
search over products of short cyclotomic integers and verifies `μ μ̄ = N`
exactly, so every "is a Hermitian norm" in §3 is proven and every "no
factor found" is only that.

### 2.4 The eigenfunction method only ever needs `K` and the `√N`

**Lemma 3.** Let `P ∈ K^{n×n}` be a Hermitian idempotent (`P^† = P = P^2`)
and `j` an index with `P_jj ≠ 0`. Then `w = P e_j` satisfies
`w_j = P_jj = ⟨w|w⟩ = N`, `N` is totally positive, and `w/√N` has the
designated component `√N`. If `w_1, …, w_{r−1}` are already chosen in the
range of `P`, the same holds for `Q = P − Σ_r |w_r⟩⟨w_r|/N_r` (the
Gram–Schmidt complement), so the multiplicity rule "take the columns in
component order, orthogonalise" never leaves `K` either.
*Proof.* `⟨w|w⟩ = ⟨e_j|P^†P e_j⟩ = ⟨e_j|P e_j⟩ = P_jj`; total positivity
because `σ(P)` is again a Hermitian idempotent for every `σ ∈ Gal(K/Q)`
(Galois conjugation commutes with complex conjugation in the abelian group
`Gal(K/Q)`), so `σ(N) = ⟨σ(w)|σ(w)⟩ > 0`. ∎

Consequences:

- **Clebsch–Gordan coefficients.** With `P = P^γ_11` on `V_α ⊗ V_β`, the
  vector `|γ 1; r⟩ ∝ w_r` and the partners `|γ k; r⟩ = P^γ_k1 w_r` have the
  same squared norm `N_r`, so **all** Clebsch–Gordan coefficients of the
  pair `(α, β)` lie in `K({√N_r})` under the real-positive rule, and in `K`
  under a norm-factor rule with `μ_r ∈ K`. Nothing beyond the finitely many
  `N_r` per group decides the field.
- **Chain-adapted vectors.** `P = Π_i P^{σ_i}_{H_i}` (commuting projectors),
  `N = P_jj`; the chain-adapted *matrices* are
  `D'(g)_ik = ⟨w_i|D(g)|w_k⟩ / √(N_i N_k)`, so their field is
  `K({√(N_i N_k) : ⟨w_i|D(g)|w_k⟩ ≠ 0 for some g})` under the real-positive
  rule and `K` under a norm-factor rule.
- **Clebsch–Gordan coefficients in a chain-adapted basis** are the
  canonical ones conjugated by `U_α ⊗ U_β` and `U_γ`, `U = [w_i/μ_i]`; their
  field is the compositum of the canonical `cg` field with the chain fields
  of the three irreps. The map's decision that the `cg` algorithm is
  basis-agnostic means the same `N`-analysis applies directly to the
  projectors written in the chain-adapted matrices.

### 2.5 The chain field must contain the subgroup characters

**Lemma 4.** If a chain-adapted basis of an irrep `γ` of `G` is adapted to
`H ⊂ G` and `σ` is a constituent of `γ ↓ H`, every field containing the
chain-adapted matrices (or vectors) contains the character values of `σ`.
*Proof.* In the chain-adapted basis `D(h)`, `h ∈ H`, is block diagonal with
blocks equivalent to `σ`, so `χ_σ(h)` is a partial trace of a matrix over
the field. ∎ Hence the conductor of a `--chain` document is at least
`lcm(M_G, cond(Q(χ_{H_1})), …)`; in practice the version-1 conductors of
the subgroups (`M_H`) are a convenient upper bound of the character
conductors, and the script works in `Q(ζ_lcm(M_G, M_{H_i}))` (the *chain
field*). Example: `O ⊃ T ⊃ C3` needs `ζ_3 ∉ Q(ζ_8)`, chain field `Q(ζ_24)`.

### 2.6 The axial families for every `n`

For `R = C_N` every irrep is one-dimensional, so every projector is `0` or
`1` on a basis vector, `N = 1`, and the Clebsch–Gordan coefficients are
`0` and `1`. For `R = D_N` (the families `Cnv`, `Dn`, `Dnd`, `Dnh`) the
irreps have dimension `≤ 2` and the canonical generator images are
`D(a) = diag(u^k, u^{−k})`, `D(b) = [[0, c], [c, 0]]` with `c ∈ {1, i}`
([Basis conventions](../spec/basis-conventions.md#canonical-basis-axial-families)).
In `E_k ⊗ E_l` the product basis vectors `e_r ⊗ e_s` are eigenvectors of
`D(a) ⊗ D(a)` with eigenvalues `u^{±k±l}`, and `D(b) ⊗ D(b)` permutes them
up to the phases `c^2 ∈ {±1}`. The projector `P^γ_11` of a two-dimensional
`γ` is a rank-one projector onto a single product vector when the four
eigenvalues are distinct (`N = 1`), and onto `(e_r⊗e_s ± c' e_{r'}⊗e_{s'})/2`
with `|c'| = 1` when two eigenvalues coincide or `γ` is one-dimensional
(`N = 1/2`); a two-dimensional times a one-dimensional irrep gives
`N = 1`. The same holds for the projective irreps `E_{k/2}` (same matrix
shapes with `c = i`) and for the type II groups (the parity sign is a
scalar). The chain-adapted vectors of a chain `D_N ⊃ D_m ⊃ C_m` or
`D_N ⊃ C_m` are the canonical basis vectors (`E_k ↓ D_m` irreducible) or
`(e_1 ± c' e_2)/√2` (`E_k ↓ D_m = A ⊕ B`), again `N ∈ {1, 1/2}`.

Therefore **[proven for all `n`]**: the Clebsch–Gordan coefficients and
chain-adapted vectors of the axial families have `N ∈ {1, 1/2}`; under the
real-positive rule the field is `K(√2) = Q(ζ_lcm(M, 8))`; `2 = |1 − i|^2`
is a Hermitian norm from `K` whenever `4 | M` (every dihedral canonical
conductor `lcm(2N, 4)`), so a norm-factor rule keeps the version-1
conductor. The survey (§3.1) confirms this on `C3v`, `D4h`, `D5h`, `D6h`,
`D7h`, `D8h`, `D9d`.

---

## 3. Computational evidence

### 3.0 The script

`docs/research/scripts/cyclotomic_closure.py` (`uv run python3
docs/research/scripts/cyclotomic_closure.py [GROUP …]`, standard library
plus `pointgroup_tables.field`, no new dependency) computes, for every
group named, the irreps of `build_model` (canonical basis, ordinary and
projective) and:

- the squared norms `N` of the Clebsch–Gordan vectors of every irrep pair
  `(α, β)` and every constituent `γ` (first partner, every multiplicity
  copy, Gram–Schmidt in component order; the designated component is
  checked to equal `⟨w|w⟩` exactly);
- for the proper polyhedral groups, the squared norms of the chain-adapted
  vectors of the chains listed in §3.3 (subgroups embedded as element sets
  of `G†` by generator words, verified to be injective homomorphisms;
  isotypic projectors in the chain field; rank of every chain projector
  checked to be `0` or `1`, their sum checked to be the dimension), and the
  products `N_i N_k` under the non-zero matrix elements `⟨w_i|D(g)|w_k⟩`;
- for every distinct `N`: whether `√N ∈ K` (`real_sqrt`: the root is an
  algebraic integer after clearing denominators, so it is an integer vector
  over the power basis; all `2^{φ(M)/2−1}` sign patterns of its real
  embeddings are tried and the candidate is verified exactly, a proof in
  both directions); otherwise whether `K(√N)/Q` is abelian (Lemma 1, with
  `real_sqrt` in `K`); if abelian, the smallest multiple `M'` of `M` with
  `√N ∈ Q(ζ_M')` (membership reduced to square roots in `K` through the
  quadratic subextensions `K(y_ψ)` of `Q(ζ_M')/K`, `y_ψ` the Gauss-type
  sums of the order-two characters of `Gal(Q(ζ_M')/K)`, searched up to
  `φ(M') ≤ 100`); and a Hermitian-norm factor `μ ∈ K` with `μ μ̄ = N`:
  `μ = c · Π ν_i` with `ν_i` short sums of roots of unity (sums of up to
  five distinct roots, sums with a doubled constant term, and `u + i v` for
  short real sums `u`, `v` when `i ∈ K`) and `c ∈ K^+`. The classes of `N`
  and of the norms `|ν|^2` modulo squares are separated by quadratic
  characters at prime ideals of residue degree one (a homomorphism killing
  squares), so `N = Π |ν_i|^2 · c^2` is a linear system over `F_2`; a
  solution is verified exactly (`c` by `real_sqrt`), and the number of
  characters is doubled while a solution fails to verify. "No factor" means
  that `N` is outside the span of these norms modulo squares, which proves
  nothing about `N` being a norm; every "factor" is an exact identity.

The whole default survey runs in about two minutes on a laptop; its output
is kept as `docs/research/cyclotomic-closure-survey.txt`.

### 3.1 Clebsch–Gordan coefficients, canonical basis

Every irrep pair `(α, β)` with `α ≤ β` in canonical irrep order, every
constituent `γ`, every multiplicity copy. "Closure" is the smallest
cyclotomic field holding every `√N` (the field of the real-positive rule);
"factor" is whether a Hermitian-norm factor `μ ∈ K` with `μ μ̄ = N` was
exhibited for every `N` (the field of a norm-factor rule is then `K`).

| Group | `M` | pairs | vectors | distinct `N` | values of `N` | closure | factor in `K` |
| --- | ---: | ---: | ---: | ---: | --- | --- | --- |
| `C3v` | 12 | 21 | 27 | 2 | `1, 1/2` | `Q(ζ_24)` | all |
| `D4h` | 8 | 105 | 142 | 2 | `1, 1/2` | `Q(ζ_8)` | all |
| `D5h` | 20 | 91 | 150 | 2 | `1, 1/2` | `Q(ζ_40)` | all |
| `D6h` | 12 | 171 | 252 | 2 | `1, 1/2` | `Q(ζ_24)` | all |
| `D7h` | 28 | 153 | 264 | 2 | `1, 1/2` | `Q(ζ_56)` | all |
| `D8h` | 16 | 253 | 394 | 2 | `1, 1/2` | `Q(ζ_16)` | all |
| `D9d` | 36 | 300 | 476 | 2 | `1, 1/2` | `Q(ζ_72)` | all |
| `T` | 12 | 28 | 44 | 4 | `1, 1/2, 1/3, 2/3` | `Q(ζ_24)` | all |
| `Td` | 8 | 36 | 74 | 8 | `1, 1/2, 1/4, 1/3, 2/3, 1/6, 3/4, 3/8` | `Q(ζ_24)` | all |
| `O` | 8 | 36 | 74 | 8 | same as `Td` | `Q(ζ_24)` | all |
| `Oh` | 8 | 136 | 272 | 8 | same as `O` | `Q(ζ_24)` | all |
| `I` | 20 | 45 | 133 | 28 | see below | `Q(ζ_840)` | all |
| `Ih` | 20 | 171 | 494 | 30 | those of `I` and `5/12`, `9/20` (pairs `(Hg, Gu)`: the pair order swaps the product components, so the designated `j` and its `N` change) | `Q(ζ_840)` | all |

The 28 squared norms of `I` (`M = 20`, `K^+ = Q(√5, √((5+√5)/2))`) are
**all rational**: `1, 1/2, 1/4, 9/16, 1/3, 2/3, 3/4, 1/9, 2/9, 1/5, 2/5,
3/5, 4/5, 9/25, 1/10, 3/10, 7/10, 9/10, 1/20, 1/6, 1/12, 1/15, 2/15, 4/15,
8/15, 1/7, 2/7, 7/20`. Their square roots need `√2 ∈ Q(ζ_8)`,
`√3 ∈ Q(ζ_12)`, `√5 ∈ Q(ζ_5)` and `√7 ∈ Q(ζ_28)`; the real-positive closure
is `Q(ζ_lcm(20, 8, 12, 28)) = Q(ζ_840)` (degree 192). The prime `7` enters
only through the pairs with multiplicity two (`I5/2 × I5/2 ⊃ T1` twice,
`T1 × I5/2 ⊃ I5/2` twice), i.e. through the multiplicity rule of the survey
(Gram–Schmidt in component order); the diagonal entry `P_jj` of a rank-two
projector is a sum of two `|u_j|^2` and need not have the prime support of
the group order. A multiplicity rule that splits `α × α` into its symmetric
and antisymmetric parts would change these `N`; this is the CG ticket's
decision, and the survey shows that it affects the field.

Observed on every group: the `N` of the type II groups (`D4h`, `D6h`,
`D8h`, `D9d`, `Oh`, `Ih`) and type III groups (`C3v`, `D5h`, `D7h`, `Td`)
are exactly those of their rotation groups, as the parity is a scalar on
every irrep.

### 3.2 Hermitian-norm factors for the Clebsch–Gordan norms

The factors the script exhibited (all verified exactly, `μ μ̄ = N`); they
are the phases a norm-factor rule would use. The axial `1/2` is
`|(1 − i)/2|^2` written in the field's own power basis.

| Field | `N` | `μ` |
| --- | --- | --- |
| `Q(ζ_8)` (`D4h`, `O`, `Td`, `Oh`) | `1/2`, `1/4` | `√2/2`, `1/2` (square roots in `K`) |
| | `1/3` | `(1 + i√2)/3` |
| | `2/3` | `(√2 + 2i)/3` |
| | `3/4` | `(1 + i√2)/2` |
| | `1/6` | `√2/6 + i/3` |
| | `3/8` | `√2/4 + i/2` |
| `Q(ζ_12)` (`T`, `C3v`, `D6h`) | `1/3` | `√3/3` |
| | `1/2` | `(ζ − ζ^2 − ζ^3)/2` (a unit-modulus multiple of `(1 − i)/2`) |
| | `2/3` | `(2 − ζ − ζ^2 − ζ^3)/3` |
| `Q(ζ_20)` (`D5h`) | `1/2` | `(1 − ζ^2 − ζ^3 + ζ^4 − ζ^6)/2` |
| `Q(ζ_28)` (`D7h`) | `1/2` | `(−ζ^4 − ζ^11)/2` |
| `Q(ζ_36)` (`D9d`) | `1/2` | `(ζ^2 − ζ^5 − ζ^8)/2` |
| `Q(ζ_16)` (`D8h`) | `1/2` | `√2/2` |

For `I` and `Ih` (`Q(ζ_20)`, `φ = (1 + √5)/2`): the rational norms factor
through `2 = |1 − i|^2`, `5 = (√5)^2`, `3 = φ^2 + φ^{−2} = |φ + i φ^{−1}|^2`
and `7 = φ^4 + φ^{−4} = |φ^2 + i φ^{−2}|^2` (the Lucas numbers `L_2`, `L_4`;
the second was found by hand after the script's first candidate families,
sums of distinct roots of unity, had missed it, which is why the `u + i v`
family of `_short_sums` exists; the script now exhibits all 28 factors of
`I` and all 30 of `Ih`). Examples from the survey: `1/3 = |(1 + √5 − i +
i√5)/6|^2`, `1/6 = |(1 + i√5)/6|^2`, `1/7 = |(3 + √5 − 3i + i√5)/14|^2`,
`2/7 = |(√5 + 3i)/7|^2`, `7/10 = |1/2 + 3i√5/10|^2`, `3/10 = |√5/10 + i/2|^2`.

The theory of §2.3 predicts these: `Q(ζ_12)`, `Q(ζ_20)`, `Q(ζ_24)`,
`Q(ζ_28)`, `Q(ζ_36)` are unramified over their real subfields, so a
totally positive rational `p/q` is a Hermitian norm iff every prime with
odd exponent in `pq` has its primes of `K^+` split in `K` or of even
valuation; `2` (`|1 − i|^2`), `3` (`|1 − ζ_3|^2` when `3 | M`, or
`φ^2 + φ^{−2}` in `Q(√5, i)`), `5` (`√5^2`) and `7` (`φ^4 + φ^{−4}` in
`Q(√5, i)`, with `φ` the golden ratio) are Hermitian norms in the fields
where they occur. For the prime-power conductor `8` the ramified prime
above `2` needs the Hilbert symbol, but the exhibited factors settle every
case that occurs.

### 3.3 Chain-adapted bases

Chains of the proper polyhedral groups, every irrep (ordinary and
projective). Subgroups: `D_n = ⟨c, c'⟩` with `c` the named axis and `c'`
the first binary rotation inverting it, `T = ⟨C3^1(1), c'⟩` with `c'` the
first binary rotation closing on order `12`, `D2` and `C2` the first binary
rotations inside the previous subgroup. Every chain projector had rank `0`
or `1` (all chains canonical). "Vectors" is the closure for the
chain-adapted vectors `w_i/√N_i`, "matrices" for the chain-adapted matrices
(`√(N_i N_k)` under the non-zero matrix elements), "factor" whether a
Hermitian-norm factor was exhibited in the chain field for every `N`
(vectors) and every `N_i N_k` (matrices).

| Chain | chain field | distinct `N` | non-cyclotomic `√N` | vectors | matrices | factor |
| --- | --- | ---: | --- | --- | --- | --- |
| `T ⊃ C3` | `Q(ζ_12)` | 1 | — | `Q(ζ_12)` | `Q(ζ_12)` | all |
| `T ⊃ D2 ⊃ C2` | `Q(ζ_12)` | 4 | `(3 ± √3)/6` (`E1/2`) | none | `Q(ζ_24)` | **vectors: 2 of 4** (`(3 ± √3)/6` not a norm from `Q(ζ_12)`, §2.3); matrices: all |
| `O ⊃ D4 ⊃ C4` | `Q(ζ_8)` | 1 | — | `Q(ζ_8)` | `Q(ζ_8)` | all |
| `O ⊃ T ⊃ C3` | `Q(ζ_24)` | 9 | `(3 ± √3)/6` (`E1/2`), `(3 ± √3)/12` (`G3/2`) | none | `Q(ζ_24)` | all |
| `O ⊃ D3 ⊃ C3` | `Q(ζ_24)` | 10 | same, and `N_i N_k = (3 ± √3)/48` (`G3/2`) | none | none | all |
| `O ⊃ T ⊃ D2 ⊃ C2` | `Q(ζ_24)` | 2 | — | `Q(ζ_24)` | `Q(ζ_24)` | all |
| `I ⊃ D5 ⊃ C5` | `Q(ζ_20)` | 2 | — | `Q(ζ_40)` | `Q(ζ_40)` | all |
| `I ⊃ T ⊃ C3` | `Q(ζ_60)` | 27 | 15 of the 27 (`1/4 − √5/30` and degree-8 elements of `Q(ζ_60)^+`; the ordinary `G`, `H` and all spinor irreps); 18 of the 59 products | none | none | all (27 vectors, 59 products) |
| `I ⊃ T ⊃ D2 ⊃ C2` | `Q(ζ_60)` | 19 | 10 of the 19 (degree-8 elements of `Q(ζ_60)^+`, spinor irreps only); 1 of the 37 products (`1/20 + √5/200`) | none | none | all (19 vectors, 37 products) |
| `I ⊃ D3 ⊃ C3` | `Q(ζ_60)` | — | **not canonical**: `H ↓ D3 = A1 ⊕ 2E`, `I5/2 ↓ D3† = … ⊕ 2E1/2`, detected by the rank check; excluded from the survey | | | |

Observations:

- A chain that ends in the cyclic group of the **principal element** of the
  canonical basis (`T ⊃ C3`, `O ⊃ D4 ⊃ C4`, `I ⊃ D5 ⊃ C5`) is already
  diagonal in the canonical basis: its chain-adapted vectors are the
  canonical basis vectors (`N = 1`) except where two eigenvalues of the
  principal element coincide, which happens only for `I5/2` of `I†`
  (`N = 1/2`, [Basis conventions](../spec/basis-conventions.md#canonical-basis-the-principle)).
- The irrational, non-cyclotomic norms come from the spinor irreps and from
  a subgroup whose axis is not the principal axis: `(3 ± √3)/6 = cos²(θ/2)`,
  `sin²(θ/2)` for `cos θ = 1/√3`, the angle between a three-fold and a
  two-fold axis of the tetrahedron. They are exactly the quantities a
  half-angle rotation of a spinor produces, and they are **provably not in
  any cyclotomic field** (Lemma 1 run exactly). The ordinary irreps of `T`
  and `O` gave `N ∈ {1, 1/2, 1/3, (2 ± √3)/6}`, all with cyclotomic square
  roots; for `I ⊃ T ⊃ C3` the ordinary `G` and `H` carry non-cyclotomic
  norms too (`G ↓ T ⊃ A`, `H ↓ T ⊃ 1E, 2E`), so the phenomenon is not
  confined to the spinor irreps.
- For every chain, every `N` was a Hermitian norm from the chain field
  enlarged to contain `√2` and `√3` (`Q(ζ_24)` for `T` and `O`; `Q(ζ_60)`
  itself for the chains of `I` through `T`), so a
  norm-factor phase rule keeps every chain-adapted vector and matrix of
  these chains cyclotomic. For `T ⊃ D2 ⊃ C2` the chain field `Q(ζ_12)`
  itself is provably insufficient (§2.3): the conductor rule must add the
  factor that makes the obstructing primes split or ramify (`8`, giving
  `√2`, does it here).
- The chain-adapted *matrices* with the real-positive rule are cyclotomic
  for `O ⊃ T ⊃ C3` and `T ⊃ D2 ⊃ C2` (the products `N_i N_k` pair the
  Galois-conjugate norms, `(3 + √3)(3 − √3) = 6`), but **not** for
  `O ⊃ D3 ⊃ C3` (`G3/2`), so the real-positive rule cannot be rescued by
  emitting only matrices.

### 3.4 What the computation did not cover

- Axial groups with `n ≥ 10`, and the `altmann-herzig` basis convention
  (its polyhedral conductors `24`, `24`, `120` make the sign search of
  `real_sqrt` cost `2^15` patterns per call; the canonical basis is what the
  map fixes as the default and the reference bases are reached from it by
  a unitary intertwiner over the convention's field).
- Clebsch–Gordan coefficients of the type II groups were computed, but they
  are the proper group's numbers (the parity is a scalar), as the tables
  show; `Td`, `C3v`, `D5h`, `D7h`, `D9d` (type III) are the rotation group's
  numbers too.
- Chains of the type II/III groups (the map embeds subgroups by rule; the
  projectors are those of the rotation-group chain with the parity
  characters attached, so the norms are the same).
- Chains of `I` other than `I ⊃ D5 ⊃ C5`, `I ⊃ T ⊃ C3`, `I ⊃ T ⊃ D2 ⊃ C2`
  (`I ⊃ D3 ⊃ C3` was found non-canonical and dropped), chains of `O` and `T`
  other than the six listed, and every chain of the axial families (covered
  by the proof of §2.6).
- Clebsch–Gordan coefficients *in* a chain-adapted basis (§2.4 reduces
  them to the compositum; not computed).
- The exact local decision "is `N` a Hermitian norm" where the search found
  no factor (none remained after enlarging the field where §3.3 says so).

---

## 4. Answers to the ticket's items

**Item 1 — which square roots occur, and are they cyclotomic?**

- *Axial families, every `n`* **[proven, §2.6]**: the only normalisation
  constant is `√2` (`N = 1/2`), in Clebsch–Gordan coefficients and in
  chain-adapted bases alike; `√2 ∈ Q(ζ_8)`, so the smallest field is
  `Q(ζ_lcm(M, 8))`: `M` itself when `8 | M` (`D4h`, `D8h`, `C8v`, …),
  `2M` when `4 ‖ M` (`C3v`, `D6h`: `24`; `D5h`: `40`; `D7h`: `56`; `D9d`:
  `72`). No non-cyclotomic radical can occur.
- *Polyhedral families, Clebsch–Gordan coefficients in the canonical
  basis* **[observed on `T`, `Td`, `O`, `Oh`, `I`, `Ih`]**: every `N` is
  rational; the radicals are `√2`, `√3`, `√6` (`T`, `O`: smallest field
  `Q(ζ_24)`) and `√2`, `√3`, `√5`, `√7` and their products (`I`: `√7` only
  through the multiplicity-two pairs under the survey's multiplicity rule;
  smallest field `Q(ζ_840)` with it, `Q(ζ_120)` without the pairs of
  multiplicity two). All cyclotomic. The nested radicals of the ticket's
  examples do not appear in Clebsch–Gordan coefficients because the
  squared norms are rational; whether that is a theorem is open
  **[conjectured]**: the diagonal entries `P_jj` of the projectors are
  `(d_γ/|G|) Σ_g χ_γ(g)^* D_α(g)_{jj} D_β(g)_{jj}` and nothing forces them
  to be rational in general, but on the surveyed groups they are.
- *Polyhedral families, chain-adapted bases* **[observed on nine canonical chains of
  `T`, `O`, `I`; the non-cyclotomicity proven exactly]**: chains ending in
  the principal element's cyclic group give `N ∈ {1, 1/2}`; chains through
  a non-principal axis give irrational `N` such as `(2 ± √3)/6` (cyclotomic
  root) and, for the spinor irreps, `(3 ± √3)/6`, `(3 ± √3)/12`, `(3 ± √3)/48`
  whose square roots **lie in no cyclotomic field** (`K(√N)/Q` is not
  Galois: `(3 + √3)(3 − √3) = 6` is not a square in `Q(√3)`). For `I` the
  chains through `T` give, already for the ordinary irreps `G` and `H` of
  `I ⊃ T ⊃ C3`, squared norms of degree eight over `Q` (elements of
  `Q(ζ_60)^+` with `√3`, `√5` and `cos(2π/15)` entangled) whose square roots
  are non-cyclotomic: `25` of the `46` distinct values across the two
  chains. This
  is basis-independent in the sense the ticket asks: the modulus of a
  chain-adapted vector's component relative to the canonical basis is a
  half-angle cosine, and no choice of canonical phases changes a modulus.

**Item 2 — does the phase freedom keep every value cyclotomic?**

- **Yes for every value met**, with the field enlarged as §3 states: for
  every squared norm of every surveyed group and chain a unit-modulus
  factor `μ/|μ|`, `μ` in the (enlarged) chain field, puts the vector into
  that field, exhibited exactly (`μ μ̄ = N`). For the chain-adapted
  matrices the same factors give entries `⟨w_i|D(g)|w_k⟩/(μ̄_i μ_k)` in the
  chain field.
- **The phase freedom is necessary**: for `T ⊃ D2 ⊃ C2`, `O ⊃ T ⊃ C3`,
  `O ⊃ D3 ⊃ C3` the real-positive rule leaves every cyclotomic field
  (vectors), and for `O ⊃ D3 ⊃ C3` even the matrices do **[proven by exact
  computation]**.
- **It is not free of charge**: `N` must be a Hermitian norm from the
  field, which the version-1 field does not always grant (`(3 ± √3)/6` is
  not a norm from `Q(ζ_12)`, §2.3, but is from `Q(ζ_24)`). The general
  statement "every `N` of every canonical chain is a Hermitian norm from
  some cyclotomic field" is **[conjectured]**; it is exactly the statement
  that the chain-adapted basis has a unitary realisation over a cyclotomic
  field, and §2.3 reduces it to a finite local computation per `N`
  (Hasse norm theorem), which the implementation can run as a check.
- For Clebsch–Gordan coefficients in the canonical basis the phase freedom
  is not needed for cyclotomicity (the real-positive closure exists) but it
  shrinks the field dramatically: `Q(ζ_20)` instead of `Q(ζ_840)` for `I`
  (the `7` of the multiplicity-two pairs is `|φ^2 + i φ^{−2}|^2`, §3.2).

**Item 3 — conductor rules for `--chain` and `cg`; do the version-1 fields suffice for `cg` in the canonical basis?**

- `cg`, canonical basis, **real-positive rule**: the version-1 fields do
  **not** suffice (`√2 ∉ Q(ζ_12)` already for `C3v`). The rule is
  `M_cg = lcm(M, 8)` for the axial families (proven for all `n`),
  `lcm(M, 24) = 24` for `T`, `Th`, `Td`, `O`, `Oh` (observed, complete for
  these groups since every pair was computed) and `lcm(M, 8, 12, 28) = 840`
  for `I`, `Ih` under the component-order multiplicity rule (`120` if the
  multiplicity rule keeps `7` out).
- `cg`, canonical basis, **norm-factor rule**: the version-1 field
  suffices on every surveyed group (`M = 8, 12, 16, 20, 28, 36` as listed),
  and §2.3 explains why (`2`, `3`, `5`, `7` are Hermitian norms from these
  fields).
- `--chain`: the conductor is at least `lcm(M_G, M_{H_1}, …)` (proven,
  §2.5); the real-positive rule is impossible in general (item 1), so the
  rule must be a norm-factor rule and the conductor must make every `N` a
  Hermitian norm: observed sufficient are `Q(ζ_24)` for every chain of `T`
  and `O` (`lcm(M_G, M_H, 24)`) and `Q(ζ_60) = Q(ζ_lcm(20, 12))` for the
  two canonical chains of `I` through `T` (no enlargement was needed there:
  `Q(ζ_60)` already contains `i`, `√3`, `√5`); `Q(ζ_120)` is the uniform
  `lcm(M_G, M_H, 24)` rule for `I`.

---

## 5. Recommendation

1. **Keep one cyclotomic field per document; do not extend the value
   model.** Nothing in the survey needs a value outside a cyclotomic field
   once the phase rule is a norm-factor rule, and §2.4 shows that the only
   candidates for leaving the field are the normalisation constants, which
   the phase rule controls. A "quadratic extension" of the value model
   (values `x/√N`) would be needed only to keep a real-positive phase rule
   for chain-adapted bases; the price (a second kind of value in every
   consumer, and `√N` that is not even cyclotomic) is not worth the
   convention, which no reference work shares anyway (Chen and Ping's own
   program writes `√(p/q) exp(i arccos(m/n))`, a unit-modulus algebraic
   phase, not a real-positive component).
2. **Phase rule (for the basis-rule ticket #43 and the CG ticket #46): a
   norm-factor rule.** Normalise `w` to `w/μ` with `μ ∈ F`, `μ μ̄ = N`,
   chosen canonically, e.g. "`μ = √N` when `√N ∈ F^+`; otherwise `μ = c·ν`
   for the first `ν` in a fixed enumeration of short sums of roots of unity
   (by length, then by exponents) with `N/|ν|^2 = c^2` a square in `F^+`,
   `c > 0`". The designated component of the normalised vector is then
   `N/μ = μ̄`, a canonical element of modulus `√N` instead of the real
   number `√N`. The rule is deterministic, reproducible from the document
   field alone, and reduces to the real-positive rule whenever `√N ∈ F`
   (every axial group with `8 | M`, `D4h`, `O ⊃ D4 ⊃ C4`, …). The reference
   comparison of the verification chapter must compare up to a diagonal
   unitary (phases), which it needs anyway for the papers' conventions.
3. **Conductor rule for `cg` (canonical basis).** With the norm-factor rule
   the version-1 conductor `M` is kept — observed sufficient on every
   surveyed group and explained by §2.3; the implementation should assert
   the Hermitian-norm factorisation for every `N` it meets (it computes
   `μ` anyway) and fall back to the closure `lcm(M, cond(√N))` with a
   `known_discrepancy`-style note if a factor is ever missing (none is
   expected: the primes of the `N` are those of `|G|` plus the ones the
   multiplicity rule introduces, and each has been checked). If the CG
   ticket prefers the real-positive rule after all, the conductor is
   `lcm(M, 8)` for the axial families (proven) and `lcm(M, 24)` for
   `T`-, `O`-type groups and `lcm(M, 840)` for `I`-type groups (observed),
   and the multiplicity rule should be chosen to avoid `7` (symmetric and
   antisymmetric split of `α × α`, a second class operator for `α ≠ β`).
4. **Conductor rule for `--chain`.** `M_chain = lcm(M_G, M_{H_1}, …, M_{H_k}, 8·3)`
   for the polyhedral groups, i.e. the chain field of §2.5 enlarged by
   `√2` and `√3` (`24` for `T` and `O` chains, `120` for `I` chains, where
   `60` was observed sufficient for the chains through `T`), and
   `lcm(M_G, M_{H_i})` for the axial families (already a multiple of `4`,
   and `N ∈ {1, 1/2}` with `2 = |1 − i|^2`). This is the smallest rule
   consistent with the survey (`T ⊃ D2 ⊃ C2` proves `lcm(M_G, M_H) = 12`
   insufficient); the implementation should again assert the norm
   factorisation of every `N` and report a missing factor as an error of
   the chain rather than silently emitting a non-cyclotomic value. The
   `chain` table and the chain-adapted `matrices` then stay in
   `Q(ζ_M_chain)`; `cg --chain` lives in the compositum, i.e. the same
   `M_chain` (the canonical `cg` field is `M | M_chain` under the
   norm-factor rule).
5. **What stays conjectural and how to retire it.** (a) "Every `N` of
   every canonical chain of every point group is a Hermitian norm from the
   recommended field" is verified on the surveyed chains and reduces, per
   `N`, to the Hasse local conditions of §2.3; the arithmetic ticket (#44)
   can implement that local test exactly (prime ideals above `2`, `3`, `5`,
   `7` of `K^+`) if an unconditional guarantee is wanted, or simply treat a
   failed factor search as a hard error. (b) "The Clebsch–Gordan squared
   norms of the polyhedral groups are rational" is observed on all six
   groups; it is not needed for the recommendation. (c) The axial statement
   is proven for all `n`.

---

## References

- R. Brauer, *On the representation of a group of order g in the field of
  the g-th roots of unity*, Amer. J. Math. 67 (1945) 461–471.
- J.-P. Serre, *Linear Representations of Finite Groups*, GTM 42, Springer
  1977, §12.3 Theorem 24 (Brauer's theorem in exponent form; statement as
  quoted in [arXiv:2409.02589](https://arxiv.org/pdf/2409.02589)).
- W. F. Reynolds, *Projective representations of finite groups in
  cyclotomic fields*, Illinois J. Math. 9 (1965) 191–198,
  <https://doi.org/10.1215/ijm/1256067877>.
- J. L. Alperin, T.-N. Kuo, *The exponent and the projective
  representations of a finite group*, Illinois J. Math. 11 (1967),
  <https://doi.org/10.1215/ijm/1256054561>.
- GroupNames (T. Dokchitser), character tables with Schur indices:
  [SL(2,3)](https://people.maths.bris.ac.uk/~matyd/GroupNames/1/SL(2,3).html),
  [CSU(2,3)](https://people.maths.bris.ac.uk/~matyd/GroupNames/1/CSU(2,3).html),
  [SL(2,5)](https://people.maths.bris.ac.uk/~matyd/GroupNames/97/SL(2,5).html),
  [Q8](https://people.maths.bris.ac.uk/~matyd/GroupNames/1/Q8.html).
- Kronecker–Weber theorem,
  [Wikipedia](https://en.wikipedia.org/wiki/Kronecker%E2%80%93Weber_theorem).
- H. Hasse, *Beweis eines Satzes und Widerlegung einer Vermutung über das
  allgemeine Normenrestsymbol*, 1931 (Hasse norm theorem;
  [Wikipedia](https://en.wikipedia.org/wiki/Hasse_norm_theorem)).
- M. Kirschmer, *Definite quadratic and hermitian forms with small class
  number*, Habilitationsschrift, RWTH Aachen 2016,
  <http://www.math.rwth-aachen.de/~Markus.Kirschmer/papers/herm.pdf>
  (Theorem 2.4.1 local–global principle for Hermitian spaces, attributed to
  Minkowski, Hasse, Landherr, Kneser, Springer; Theorem 3.2.3 local
  classification by rank and determinant modulo norms; Proposition 2.1.5).
- W. Scharlau, *Quadratic and Hermitian Forms*, Grundlehren 270, Springer
  1985, Chapter 10 (the proof of the local–global principle Kirschmer
  refers to).
- L. C. Washington, *Introduction to Cyclotomic Fields*, 2nd ed., GTM 83,
  Springer 1997, Proposition 2.8 (`1 − ζ_n` is a unit for `n` not a prime
  power; as cited in [arXiv:2407.02002](https://arxiv.org/html/2407.02002)).
- J. M. Masley, H. L. Montgomery, *Cyclotomic fields with unique
  factorization*, J. Reine Angew. Math. 286/287 (1976) 248–256 (class
  number one for `M = 3, 4, 5, 7, 8, 9, 11, 12, 13, 15, 16, 17, 19, 20, 21,
  24, 25, 27, 28, 32, 33, 35, 36, 40, 44, 45, 48, 60, 84`; list as quoted in
  [arXiv:1507.07520](https://arxiv.org/pdf/1507.07520), Theorem 72). Not
  used in a proof above; it is why a norm factor `μ` can be taken integral
  up to a rational denominator in every document field but `Q(ζ_120)`.
- K. Hymabaccus, D. Pasechnik, *RepnDecomp: A GAP package for decomposing
  linear representations of finite groups*, J. Open Source Softw. 5(50)
  (2020) 1835, <https://arxiv.org/abs/2007.02459>; manual chapter 3
  <https://docs.gap-system.org/pkg/repndecomp/doc/chap3_mj.html>
  (`UnitaryRepresentation` via `LDLDecomposition`).
- GAP Forum thread *Repsn: constructing representations with real
  coefficients* (2021),
  <https://www.gap-system.org/ForumArchive2/2021/006241.html>
  (Repsn does not compute unitary representations; RepnDecomp does, then
  "finding an invariant quadratic form and diagonalising").
- P. van den Broek, J. F. Cornwell, *Clebsch–Gordan coefficients of symmetry
  groups*, phys. stat. sol. (b) 90 (1978) 211–224,
  <https://doi.org/10.1002/pssb.2220900123> (Clebsch–Gordan vectors as
  eigenvectors of projection matrices, orthonormalised).
- R. Dirl, *Clebsch–Gordan coefficients: General theory*, J. Math. Phys. 20
  (1979) 659, <https://doi.org/10.1063/1.524128>.
- P. H. Butler, *Point Group Symmetry Applications: Methods and Tables*,
  Plenum 1981, <https://doi.org/10.1007/978-1-4613-3141-4> (jm factors and
  3j, 6j symbols with a chapter *Notation and Phase Conventions*; the
  `I ⊃ D5 ⊃ C5` tables Fan, Chen and Draayer compare with).
- J.-Q. Chen, M.-J. Gao, G.-Q. Ma, *The representation group and its
  application to space groups*, Rev. Mod. Phys. 57 (1985) 211 (the
  eigenfunction method).
- J.-Q. Chen, P.-D. Fan, *Algebraic solutions for point groups: The
  tetrahedral group for the group chain T ⊃ D2 ⊃ C2*, J. Math. Phys. 39
  (1998) 5502–5518, <https://doi.org/10.1063/1.532547>.
- P.-D. Fan, J.-Q. Chen, J. P. Draayer, *Algebraic expressions for
  symmetry-adapted functions of the icosahedral group in spinor space*,
  Acta Cryst. A55 (1999) 1049–1060, <https://doi.org/10.1107/S0108767399010272>
  (chain `I ⊃ D5 ⊃ C5`, double-induced technique).
- J.-Q. Chen, J.-L. Ping, *Program packages for point groups and space
  groups with subgroup chain symmetry adaptation*, Comput. Phys. Commun.
  120 (1999) 71–93, <https://doi.org/10.1016/S0010-4655(99)00201-5>
  (abstract: characters, irreducible matrices and Clebsch–Gordan
  coefficients of the crystallographic point groups adapted to canonical
  chains `G ⊃ H ⊃ K`, "with exact values in the form `√(p/q) exp(iπm/n)`
  or `√(p/q) exp(i arccos(m/n))`", i.e. a square root of a rational times a
  unit-modulus algebraic number that need not be a root of unity: the
  norm-factor shape of §2.3).
- P. W. Fowler, A. Ceulemans, *Symmetry relations in the property surfaces
  of icosahedral molecules*, Mol. Phys. 54 (1985) 767,
  <https://doi.org/10.1080/00268978500103151> (Clebsch–Gordan tables of `I`
  by the projection-matrix method).
