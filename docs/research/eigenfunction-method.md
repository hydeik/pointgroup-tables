# The eigenfunction method in the Chen–Fan–Ping papers

Research note for wayfinder ticket
[Survey of the eigenfunction method in the Chen–Fan–Ping papers](https://github.com/hydeik/pointgroup-tables/issues/39),
part of [Map: version 2, chain-adapted bases and Clebsch–Gordan coefficients](https://github.com/hydeik/pointgroup-tables/issues/38).
Vocabulary follows `GLOSSARY.md` (group chain, canonical chain, chain-adapted
basis, chain label, Clebsch–Gordan coefficient, coupling multiplicity, section,
lifted/barred element, factor system, canonical basis, principal element).
The version-1 conventions compared against in section 5 are those of
`docs/spec/double-groups.md`, `docs/spec/irrep-labels.md` and
`docs/spec/basis-conventions.md`.

Sources are the fifteen papers under `PDF/` (the user's copies, gitignored).
They are cited by the keys below; page numbers are the journal pages printed
on the PDF, equation and table numbers are the papers' own. Text was extracted
with `pypdf`; equations and tables extract badly, so every matrix quoted here
was re-checked against the surrounding prose, and places where a value could
not be read reliably are marked *(extraction uncertain)*.

| Key | Paper |
| --- | ----- |
| [P89] | J.-L. Ping, Q.-R. Zheng, B.-Q. Chen, J.-Q. Chen, *Computer generated subgroup-symmetry adapted irreducible representations and CG coefficients of space groups by the eigenfunction method*, Comput. Phys. Commun. 52 (1989) 355–373 |
| [CP97] | J.-Q. Chen, J.-L. Ping, *Algebraic expressions of irreducible bases for the molecule C60*, J. Mol. Spectrosc. 186 (1997) 403–415 |
| [CF98a] | J.-Q. Chen, P.-D. Fan, *Algebraic solutions for point groups: the tetrahedral group for the group chain T ⊃ C3*, J. Math. Phys. 39 (1998) 467–488 |
| [CF98b] | J.-Q. Chen, P.-D. Fan, *Algebraic solutions for point groups: the tetrahedral group for the group chain T ⊃ D2 ⊃ C2*, J. Math. Phys. 39 (1998) 5502–5518 |
| [CF98c] | J.-Q. Chen, P.-D. Fan, *Algebraic expressions of the Clebsch–Gordon coefficients of the point group T†*, J. Math. Phys. 39 (1998) 5519–5535 |
| [CP99] | J.-Q. Chen, J.-L. Ping, *Program packages for point groups and space groups with subgroup chain symmetry adaptation*, Comput. Phys. Commun. 120 (1999) 71–93 |
| [FC99] | P.-D. Fan, J.-Q. Chen, *Algebraic solutions for the octahedral group: group chain O ⊃ C4*, Int. J. Quantum Chem. 74 (1999) 7–22 |
| [FCMB99] | P.-D. Fan, J.-Q. Chen, L. McAven, P. Butler, *Unique Euler angles and self-consistent multiplication tables for double point groups*, Int. J. Quantum Chem. 75 (1999) 1–9 |
| [PC99] | J.-L. Ping, J.-Q. Chen, *Subgroup-chain symmetry-adapted irreducible matrices and CG coefficients of point groups*, Int. J. Quantum Chem. 75 (1999) 67–80 |
| [FCD99a] | P.-D. Fan, J.-Q. Chen, J. P. Draayer, *Reduced projection operators and algebraic expressions for the symmetry-adapted functions of the icosahedral group*, Acta Cryst. A55 (1999) 871–883 |
| [FCD99b] | P.-D. Fan, J.-Q. Chen, J. P. Draayer, *Algebraic expressions for symmetry-adapted functions of the icosahedral group in spinor space*, Acta Cryst. A55 (1999) 1049–1060 |
| [CFP00] | J.-Q. Chen, P.-D. Fan, J. Paldus, *Algebraic solutions for point groups: cubic groups G in the group chain G ⊃ T ⊃ D2 ⊃ C2*, Int. J. Quantum Chem. 76 (2000) 585–599 |
| [CFMB00] | J.-Q. Chen, P.-D. Fan, L. McAven, P. Butler, *Algebraic solutions for all dihedral groups*, J. Math. Phys. 41 (2000) 8196–8222 |
| [F01] | P.-D. Fan, *Algebraic solutions for point groups: the octahedral group for the group chain O ⊃ T ⊃ C3*, Int. J. Quantum Chem. 83 (2001) 259–270 |
| [F03] | P.-D. Fan, *Algebraic expressions of Clebsch–Gordon coefficients for the group chain O† ⊃ C4†*, Int. J. Quantum Chem. 91 (2003) 541–550 |

Notation of the papers kept in this note: `C` is the CSCO of `G`, `C(s)`
that of the subgroup chain, a bar marks an intrinsic operator, `θ = +1`
(single-valued, "s.v.") or `−1` (double-valued, "d.v.") is the representation
type index, `ξ = 0` or `1/2` the same index as an exponent, `ρ_m =
exp(−2πim/n)` the eigenvalue of `C_nz = exp(−2πi J_z/n)` on the quantum
number `m`, `ε = exp(−2πi/3)`, `|λ|` or `h_λ` the dimension of irrep `λ`.
Our own symbols: `s(g)` lifted element, `Ē` full turn, `ω` factor system,
`ζ_M = exp(2πi/M)`, `e` the exponent of a canonical-basis eigenvalue
`ζ_N^e`.

---

## 1. The method

### 1.1 Representation group and the double group

All papers treat single- and double-valued irreps at once through the
**representation group** (rep group) `𝒢`: the set of representation operators
`{R_i, −R_i}` of the double group `G† = {R_i, R̃_i}` in a spinor space, in
which `R̃_i = −R_i` is linearly dependent on `R_i`. The rep group has `|G|`
linearly independent elements and its class operators are linearly dependent
(some vanish, some differ by sign); the number of inequivalent d.v. irreps
equals the number `N` of linearly independent class operators (generalised
Burnside theorem) ([P89] §2 eqs 2.1–2.10; [CP99] §2.1 Remarks 1–2 and
Definitions 1–2; [FCD99b] §2 Definitions 1–3, Theorem 1, Examples 1–3). A
class operator of `G†` whose class mixes tilded and untilded elements is the
null operator in the rep group ([FCD99b] Example 2: `C_3 = C_4 = C_5 = 0` for
the quaternion group `Q ≅ D2†`). Multiplication in the rep group is
`R_i R_j = η(i, j) R_k` with `|η| = 1` ([P89] eq 2.6), i.e. exactly our factor
system `ω(g, h) = ±1` once the untilded set is fixed (a section).

The double-group multiplication table is built from Euler angles with the
third angle extended to `0 ≤ γ < 4π` (`0 ≤ α < 2π`, `0 ≤ β ≤ π`, `α = 0` when
`β = 0` or `π`), the matrix `D^{1/2}(α, β, γ)` of eq (11) in [FCMB99] (Rose's
convention, `α` and `γ` opposite to Bradley–Cracknell) being identified with
`cos(ω/2) − i n·σ sin(ω/2)` of the rotation by `ω ∈ [0, 4π)` about `n(θ, φ)`;
the **unique assignment** `R_i = D^{1/2}(α_i, β_i, γ_i)`, `R̃_i =
−D^{1/2}(α_i, β_i, γ_i)` removes the ambiguity of Bradley & Cracknell p. 424
([CF98a] §III eq 10; [CP99] §2.2; [FC99] §2; [FCMB99] §3 eqs 8–14). In the
resulting tables "elements `R, S, T, …` belong to one class and `R̃, S̃, T̃, …`
to a different class instead of a mix" ([CF98a] Remark 1 p. 473; [FCMB99]
Remark 1 p. 5). This is a *section* in the glossary sense; its relation to
the version-1 canonical section is in section 5.1.

### 1.2 CSCO-I, -II, -III

Let `C_i` (`i = 1..N`) be the linearly independent class operators of `𝒢`.
Because the projection operator `P^{(Γ)} = (|Γ|/|G|) Σ_i χ_i^{Γ*} C_i` is an
eigenvector of every class operator (Schur), a linear combination

```
C = Σ_i a_i C_i          (CSCO-I, or "the CSCO of G")
```

with `N` distinct eigenvalues `λ^Γ = Σ_i a_i λ_i^Γ` labels the irreps
uniquely; the eigenvalue of a single class operator is `λ_i^Γ = g_i χ_i^Γ /
h_Γ` (`g_i` class size) ([CP99] §3 eqs 7–11; [CF98a] eq 3; [PC99] eq 2). In
the class space `C` is represented by `C^{(c.s.)} = Σ_i a_i D(C_i)` with
`D_{kj}(C_i) = c_{ij}^k` the class structure constants, and
`C^{(c.s.)} χ^{(λ)*} = λ χ^{(λ)*}`: the complex conjugate of the primitive
character vector is an eigenvector of the CSCO-I ([CP99] eqs 12–14; [P89]
§2 end). Phase of the character eigenvector: `χ_1^λ = |λ| > 0`, then
normalised ([CP99] §3 p. 77).

For a canonical chain `G ⊃ H ⊃ K` with CSCOs `C`, `C_H`, `C_K` and
eigenvalues `λ, μ, ρ`, the set `(C, C_H, C_K)` is the **CSCO-II**, complete in
any multiplicity-free representation space; a chain-adapted basis vector is a
simultaneous eigenvector with eigenvalues `(λ, μ, ρ)`, and `(μ, ρ)` is the
component index of `λ` ([CP99] §4 eqs 15–16; [PC99] eqs 11–13; [CF98a] eq 6).
For computation the chain operators are combined into one,

```
C(s) = a C_H + b C_K,    m = a μ + b ρ,
```

with `a, b` chosen so that `m ↔ (μ, ρ)` is one-to-one ([CP99] eq 17; [PC99]
eqs 3–4). Adaptation to a single subgroup is the case `C(s) = C_K`.

The **intrinsic group** `Ḡ = {R̄_i}` acts in the group space by
`R̄_i R_j = R_j R_i` (multiplication from the right); it is anti-isomorphic to
`G`, commutes with it, has the same CSCO `C̄ = C`, and `{R̄_i^{−1}}` is the
right regular representation ([CP99] §5 eqs 18a–b; [CF98a] eq 5; [PC99] eqs
8–9). The **CSCO-III** `(C, C(s), C̄(s))` has `|G|` distinct sets of
eigenvalues `(λ, m, m̄)` in the group space, and its simultaneous
eigenvectors are the (complex conjugates of the) irreducible matrix element
functions: with the normalised irreducible symmetry operator (ISO)

```
P^{(λ)m̄}_m = sqrt(h_λ/|G|) Σ_a D^{(λ)}_{m m̄}(R_a)* R_a,        (CF98a eq 4; CP99 eq 19)
(C, C(s), C̄(s)) P^{(λ)m̄}_m = (λ, m, m̄) P^{(λ)m̄}_m,              (CF98a eq 6; CP99 eq 22)
```

`m` is the external and `m̄` the intrinsic quantum number ("the analogue of
`(J², J_z, J̄_z)`"); `P` is a `G ⊃ G(s)` basis vector `|λ m⟩` and a
`Ḡ ⊃ Ḡ(s)` basis vector `|λ m̄⟩` ([CF98a] p. 470). In matrix form the vector
`d^{(λ)*}_{m m̄} = (D_{m m̄}(R_1), …, D_{m m̄}(R_g))^{*}` satisfies
`(C, C(s), C̄(s)) d = (λ, m, m̄) d` where the representative matrices in the
group space are read off the group table: if the `i`-th row of the table is
`(a_1 … a_g)` then `D(R_i) = {a_1 … a_g}` (entry `±1` at row `|a_j|`, column
`j`, sign of `a_j`), and if the `k`-th column is `(b_1 … b_g)` then
`D(R̄_k) = {b_1 … b_g}` ([CP99] eqs 27–29 with the `D3†` example). Again one
operator is formed,

```
K = k_1 C + k_2 C(s) + k_3 C̄(s),    eigenvalue k = k_1 λ + k_2 m + k_3 m̄,
```

and `K` must have `|G|` distinct eigenvalues; three preset coefficient sets
are tried, otherwise the user is prompted ([CP99] eqs 25–26; [P89] eq 4.1).
Deleting `C̄(s)` gives the CSCO-II operator `Λ = k_1 C + k_2 C(s)` used for
Clebsch–Gordan coefficients ([CP99] eq 31; [P89] eq 4.4).

Concrete CSCO-I used in the papers (class sums; `θ = +1` s.v., `−1` d.v.):

| Group | CSCO-I | Eigenvalues | Source |
| ----- | ------ | ----------- | ------ |
| `T` | `C = C31⁺ + C32⁺ + C33⁺ + C34⁺` (class `4C3⁺`; not Hermitian) | `4, 4ε, 4ε*, 0` (A0, A1, A−1, F); `2, 2ε, 2ε*` (E1/2, E3/2, E−3/2) | [CF98a] eq 13, Table I; [CF98c] Table I |
| `O` | `C = Σ (C4x⁺ + C4y⁺ + C4z⁺ + C4x⁻ + C4y⁻ + C4z⁻)` (class `6C4`); the class `6C2` used in the boson paper is valid for s.v. only | `6, −6, 0, 2, −2` (A1, A2, E, F1, F2); `3√2, −3√2, 0` (E1/2, E3/2≡Γ7, G3/2) | [FC99] §3 and Table III; [CFP00] eq 1a, Table I |
| `I` | `C = Σ_{j=1}^{6} (C5,j + C5,j⁻¹)` (class `12C5`) | `12, −8c₂ = 2+2√5, −8c₄ = 2−2√5, −3, 0` (A, F1, F2, G, H), `c_k = cos(2πk/5)`; d.v. `12c_{1/2} = 3+3√5, 12c_{3/2} = 3−3√5, 3, −2` (E1/2, E3/2≡E7/2, G3/2, I5/2) | [FCD99a] eqs 3–4, Table 3; [FCD99b] eqs 12–13, Table 1 |
| `Td` | class `6S4` (`CSCO-I: 1*C(3)`, `C(3) = Σ S4^±`) | `3√2, −3√2, 0` (printed `4.243, −4.243, 0`) for Γ6, Γ7, Γ8 | [CP99] test run 1 p. 89 |
| `D4` | `((1+θ)/2)(C2x + C2y) + (2+θ)(C4z⁺ + C4z⁻)` = `C2x + C2y + 3(C4z⁺+C4z⁻)` (s.v.), `C4z⁺ + C4z⁻` (d.v.) | `8, 4, −4, −8, 0`; `√2, −√2` | [PC99] §3 p. 69, Table XII–XIII |
| `D2`, `D2'`, `C2v`, `C2v'` | `(1+θ) C2y + C2z` (and the primed variants `(1+θ)C2e + C2y`, `(1+θ)σ_y + C2z`, `(1+θ)σ_db + C2z`) | `3, 1, −1, −3` (A, B1, B2, B3), `0` (E1/2) | [PC99] Table V, VI |
| `D3`, `C3v` | `C2b + C2f + θ C2e`; `σ_db + σ_df + θ σ_de` | `3, −3, 0; 0, 3i, −3i` | [PC99] Table V, XV |
| `D6` | `5(C6z⁺ + C6z⁻) + (C21'' + C22'' + C23'')` | — | [PC99] Table V |
| `D2d`, `C4v` | as `D4` with `S4z^±`, `σ` in place of `C4z^±`, `C2` | — | [PC99] Table V |
| `Cn`, `S_n` | `C_nz`, `S_nz⁻` | `ρ_m` | [PC99] p. 69 |
| `G × Ci` | `C_G + 50 Î` (coefficient 50 separates the `±` sets) | `λ_i ± 50` | [PC99] p. 69 |

The CSCO of `T`, `O`, `I` and the cyclic groups is a single class operator;
the dihedral groups need a combination that depends on `θ` ([CF98a] p. 468;
[PC99] p. 69). The program [CP99] finds a CSCO-I automatically: diagonalise one
class operator; if it has `N` distinct eigenvalues it is the CSCO, else combine
two or three with preset coefficients (three sets) ([CP99] §3 p. 77; [P89]
§4.1(c)). Chain operators `C(s)` for every preset chain of `T, Td, O, Oh`
are tabulated with element indices ([CP99] Table 6; [PC99] Table XVI), e.g.
`T ⊃ C3`: `2C31⁺ + C32⁺ + C33⁺ + C34⁺` (`= C_T + C31⁺`); `T ⊃ C2x` (d.v.):
`C2x + Σ C3i⁺`; `T ⊃ D2` (s.v.): `2C2y + C2z + Σ C3i⁺`; `D4 ⊃ D2` (s.v.):
`2C2x + 4C2y + C2z + 6(C4z⁺ + C4z⁻)`; `D4 ⊃ C2x` (d.v.): `C2x + C4z⁺ +
C4z⁻`; `D3 ⊃ C31`: `C31⁺ + C2b + C2f + θ C2e`.

### 1.3 Clebsch–Gordan coefficients by the CSCO-II

CG coefficients are the expansion coefficients

```
ψ^{(λ'')τ}_{a''} = Σ_{a a'} C^{λ a, λ' a'}_{λ'' a'' τ} ψ^{(λ)}_a ψ^{(λ')}_{a'},   τ = 1 … (λ λ' λ''),
```

([CP99] eq 34; [P89] eq 3.2). The CSCO-II `Λ = k_1 C + k_2 C(s)` is
represented in the Kronecker product space by `D^{(λ×λ')}(Λ)` with
`D^{(λ×λ')}(R_i) = D^{(λ)}(R_i) ⊗ D^{(λ')}(R_i)`; its eigenvalue
`ℓ = k_1 λ'' + k_2 a''` identifies the coupled irrep and component, the
degeneracy `d` of the eigenvalue is the coupling multiplicity `(λ λ' λ'')`,
and the `d` orthogonal eigenvectors `U^{(λ'')τ}_{a''}` (`τ = 1 … d`) are the
CG vectors for one component `a''` ([CP99] eqs 35–39; [P89] eqs 3.5–3.7). The
other components are **not** computed by further diagonalisation but by the
irrep matrices, so that CG columns of one `λ''` have the same relative phases
as the irrep matrices: choose an element `R_i` whose column `a''` of
`D^{(λ'')}` has exactly one non-vanishing off-diagonal entry `D_{b'' a''}(R_i)`;
then

```
U_{b''} = [D^{(λ'')}_{b'' a''}(R_i)]^{−1} [D^{(λ×λ')}(R_i) − D^{(λ'')}_{a'' a''}(R_i)] U_{a''};
```

if no such element exists, two elements with two off-diagonal entries and a
2×2 linear system ([CP99] eqs 40–42; [P89] eqs 3.8–3.10). The three irreps may
be adapted to different chains; the code assumes the same chain for the
factors in the s.v.×s.v. and d.v.×d.v. cases and for the two d.v. irreps in
the s.v.×d.v. case ([CP99] §8.3). Multiplicity separation is "equally simple
disregard what type it belongs to" ([P89] §5 p. 365) — i.e. the degenerate
eigenspace is handed to the eigensolver, whose orthonormal basis is arbitrary;
the papers give no rule for `τ`.

### 1.4 The double-induced technique (algebraic solutions for `G ⊃ C_n`)

For a chain whose subgroup is cyclic, `G(s) = C_n = ⟨C_nz⟩`, [CF98a] §II
replaces the `|G|`-dimensional eigenproblem by one of dimension `q`, the
number of double cosets. Steps (quoted from [CF98a] pp. 470–471):

1. Decompose `G = Σ_{i=1}^{q} C_n k_i C_n` into double cosets with fixed
   generators `k_i` (`k_1 = e`); the first `p` coincide with one-sided cosets
   (`k_i C_n = C_n k_i`).
2. Form the **double-induced representation** by sandwiching the generators
   between projection operators of `C_n` and of the intrinsic `C̄_n`:
   `φ^i_{ρρ̄} = const · P_ρ k_i P_ρ̄`, `P_ρ = Σ_k (ρ*)^k C_nz^k`. For
   `i ≤ p`, `k_i P_ρ̄ = P_{f_i(ρ̄)} k_i` with `f_i(ρ̄) = ρ̄` or `ρ̄*`, so only
   `pn + (q−p)n² = |G|` vectors are linearly independent (eq 9d); the
   double-induced representation is induced from the irrep `(ρ, ρ̄)` of
   `C_n × C̄_n` and generalises the usual (left-)induced representation
   ([CF98a] §V.A eq 51).
3. Build the `q × q` matrix `M_{ij} = ⟨φ^i|C|φ^j⟩` of the CSCO-I; only the
   coefficients of one chosen representative element per basis vector are
   needed (the `C` rows of the group table are scanned for them) ([CF98a]
   p. 476; [CF98b] p. 5512 corrects three misprints: `M̃` should read `M` in
   eqs 9g, 24, 27 of [CF98a]).
4. Diagonalise `M(ρ, ρ̄)` with `ρ, ρ̄` as parameters: `P^{(λ)ρ̄}_ρ = Σ_i u_i
   φ^i_{ρρ̄}`, `M u = λ u`. Single- and double-valued irreps come out together,
   distinguished only by the admissible values of `(λ, ρ, ρ̄)`.
5. Read the irreducible matrices of the generators off the coefficients: the
   coefficient of `P_ρ k_i P_ρ̄` in the normalised ISO is `sqrt(h_λ/|G|)
   D^{(λ)}_{ρρ̄}(k_i)*` ([CF98a] eq 8, step 5). All other elements follow from
   `D^{(λ)}_{m m̄}(C_nz^j k_i C_nz^k) = ρ_m^j D^{(λ)}_{m m̄}(k_i) ρ_m̄^k`
   ([CF98a] eq 39; [FC99] eq 34; [FCD99a] eq 26); each paper prints the `(j, k)`
   pair of every element.
6. Keep one unnormalised column `ρ̄` per irrep as the **generalized projection
   operator** `P^{(Γ)}_ρ` ([CF98a] eq 40).

The **reduced (effective) projection operator** ([CF98c] §II; [FC99] eq 27;
[FCD99a] eq 14; [F03] eqs 4–5)

```
℘^{(λ)μ̄}_μ = Σ_{i=1}^{q} M_i^{−1} d^{(λ)}_{μμ̄}(β_i)* β_i,
P^{(λ)μ̄}_μ = sqrt(h_λ/|G|) P_μ ℘^{(λ)μ̄}_μ P_μ̄,
```

with `M_i = |C_n β_i ∩ β_i C_n|` the number of times an element occurs in the
double coset (`n` for the one-sided cosets, `1` otherwise; `M_i = 3, 1`
for `T`, `4, 4, 1` for `O`, `5, 5, 1, 1` for `I`), may replace the ISO when
acting on a trial state provided the **magnetic quantum number conservation
rules** `m ≐ μ` (final state) and `m̄ ≐ μ̄` (trial state), `≐` meaning
congruence mod `n`, are imposed; `q = 2, 2, 3, 4` for `D_n, T, O, I`, so the
ISO of `I_h` has four terms ([CF98c] Remarks 1–4). "Only the irreducible
matrices of the coset generators are needed for constructing the SAFs and CG
coefficients … listing the matrices for all elements is a waste" ([CF98c]
Remark 3).

The **factorization lemma** ([CF98b] §II) covers chains that are not
`G ⊃ C_n`: for `G = H ∧ K` with `H, K` Abelian and projection operators
`P_ρ` (of `H`), `℘_ν` (of `K`), the generalized ISO `P^{(λ)}_{ρ,ν} = P_ρ ℘_ν`
is already an eigenvector of `(C(s), C̄(s'))`, and `C` is redundant (Dirac's
theorem); the conventional ISO is `P^{(λ)}_{ρρ'} = P_ρ ℘_ν P_{ρ'}` (eqs
7a–d). For one-dimensional irreps the product form holds for any semidirect
product (Remark 4). This gives `T = D2 ∧ C3` ⊃ `D2` directly ([CF98b] §III)
and all dihedral groups `D_n = C_n ∧ C2x` without any group table
([CFMB00] §II.B: "in the algebraic approach no group table or factor system is
required", §VII.B).

The **right-induced technique** ([CFP00] §2; [F01] §3) builds `O` from `T`:
decompose `O = T(e + k̂_2)` into right cosets (`k̂_2 = C4z⁻` for `O ⊃ T ⊃
D2 ⊃ C2`, `C2b` for `O ⊃ T ⊃ C3`). The vectors `φ = P^{(γ)}` (ISO of `T`) and
`ψ = P^{(γ)} k̂_2` carry the quantum numbers `(γ, ρ; γ, ρ̄)` and
`(γ, ρ; γ*, f(ρ̄))` of `(C_T, C(s); C̄_T, C̄(s))` because `k̂_2 C_T =
C_T^† k̂_2` ([CFP00] eqs 6–10). For an irrep of `T` with complex `γ` (`A±1`,
`E±3/2`) the two vectors are already irreducible bases of `O` (irreps `E`,
`G3/2`); for real `γ` (`A0`, `F`, `E1/2`) the 2×2 matrix of `C_O` in
`(φ, ψ)` is diagonalised: `M = −2ν̄ [[0,1],[1,0]]` for `F` (eigenvalues
`±2` → `F1`, `F2`), `M = 3√2 [[0, i^{ν̄}],[i^{−ν̄}, 0]]` for `E1/2`
(eigenvalues `±3√2` → `E1/2`, `E3/2` of `O`) ([CFP00] eqs 25, 41; [F01]
eqs 24, 49). Result: `P^{(F_σ)} = (1/√2)[P^{(F)} + (−1)^σ P^{(F)}_{swapped} k̂_2]`
etc. ([CFP00] eq 27, 43; [F01] eqs 26, 51). This is how "the subduction rules
are directly built into the algebraic expressions" and how the SAFs of `O`
are shown to be (symmetrised) SAFs of `T` ([F01] eqs 34, 41, 58).

### 1.5 What the eigenproblem leaves open, and how each paper fixes it

**Irrep identification and order.** The eigenvalue `λ` only labels irreps; the
conventional name is attached afterwards through the characters (`λ^Γ = g_i
χ_i^Γ/h_Γ`) ([CF98a] p. 468 "using the known characters is only for
identifying irreps"). The programs print irreps "in decreasing order of the
eigenvalues of the CSCO-I" ([P89] §4.1(d)) and then switch to Koster labels
([CP99] §3). The algebraic papers use the eigenvalue itself as the irrep
label (modified Mulliken notation, section 5.3).

**Component order.** Components are the eigenvalues of `C(s)`, printed as
Koster labels of the subgroups `[Γ_j γ_k, …]` ([CP99] §8; [PC99] p. 77), or
as the quantum number `m` of the cyclic subgroup. In the algebraic papers
rows/columns are ordered by `m` descending from `+|m|` (e.g. `E_{1/2}`:
`m = 1/2, −1/2` ([CF98b] eq 38); `G3/2` of `O ⊃ C4`: `μ = −3/2, −1/2, 1/2,
3/2` ([FC99] eq 33b); `I5/2`: `5/2, 3/2, 1/2, −1/2, −3/2, −5/2` ([FCD99b]
p. 1056); `G3/2` of `O ⊃ T ⊃ C3`: `(E3/2 3/2), (E3/2 1/2), (E−3/2 −1/2),
(E−3/2 3/2)` ([F01] eq 53); `G3/2` of `O ⊃ T ⊃ C2`: `(E3/2 1/2), (E3/2 −1/2),
(E−3/2 1/2), (E−3/2 −1/2)` ([CFP00] eq 46)). The `F` of `T ⊃ D2` is ordered
`(νν') = (+−, −+, −−) ≡ (x, y, z)` ([CF98b] eq 15c; [CFP00] Table III).

**Degenerate eigenvalues.** By construction the CSCO-I/-II/-III have distinct
eigenvalues (coefficients are re-chosen until they do). Genuine degeneracy
remains in two places: (i) coupling multiplicity (section 1.3; eigensolver
basis in the programs, trial states in the algebraic papers, section 3);
(ii) `I5/2` of `I ⊃ C5`, where `μ = ±5/2` have the same `C5` eigenvalue
`−1`: `I ⊃ C5` is not canonical for d.v. and `D5` must be inserted; the two
vectors are separated by diagonalising one binary rotation `C2,8` and its
intrinsic partner `C̄2,8` (eigenvalues `∓i` ↔ `D5` irreps `[5/2±]`, Koster
`Γ7, Γ8`), `ψ_{[5/2±]} = (1/√2)(ψ_{5/2} ± i ψ_{−5/2})` ([FCD99b] §3
Modification 3, eqs 34–37, 44).

**Phases of the eigenvectors.** The programs ([CP99] §5, after eq 29b):
(1) `D^{(λ)}_{mm}(e) = 1` fixes the diagonal functions; (2) the phase of the
functions `d_{m m̄}` with `m̄ ≠ m` for one row "can be chosen freely, or in
accordance with a convention if such a convention exists" — one free phase
per basis vector beyond the first; (3) all remaining `d_{m m̄}` are fixed by
consistency with the group table: compute `D_{m m̄}(R_i)` from the already
fixed functions through `D(R_i R_a) = D(R_i) D(R_a)` (eq 30) and demand that
the eigenvector agree. In [PC99] (abstract, Discussion p. 79) the free phases
and component orders of the cubic groups are chosen so that
`D^{(Γ_i)}(G) = D^{(Γ_{i+1})}(G) ⊗ D^{(Γ_2)}(G)` for `i = 4, 6` (`F2 = F1 ⊗
A2`, `E3/2 = E1/2 ⊗ A2`) for `G = Td, O, Oh`, and `D^{(Γ_7)}(T) =
D^{(Γ_6)}(T)*`; the d.v. `G3/2` matrices of `O ⊃ T ⊃ C3` and `O ⊃ T ⊃ C2` are
block-diagonal `E3/2 ⊕ E−3/2` on `T` and anti-block-diagonal off `T`.

The algebraic papers fix phases case by case:
- `T ⊃ C3`: the d.v. off-diagonal ISO carries `−i`, "chosen for
  convenience", so that `(P^{m̄}_m)† = P^{m}_{m̄}` and `E−3/2 = (E3/2)*`
  ([CF98a] eqs 35–36).
- `T ⊃ C2`: `√i` for `m = −1/2` so that the SAFs are real, `(√i)³` for
  `m = +1/2` from hermiticity ([CF98b] eq 35c); re-phased by `i` in [CF98c]
  eq 32 so that `D(C31^±)` are symmetric.
- `O ⊃ C4`: "phase factors determined by the requirement that the SAFs of
  `E, F0, F2` in `l = 1, 2` and `E1/2, E3/2, G3/2` in `j = 1/2, 3/2, 5/2` …
  are identical with Butler's" ([FC99] p. 13); explicit phase tables
  `θ^1, θ^2 = −1` for `(μ, μ̄) = (1/2, −1/2), (3/2, −3/2)` and `θ^3 = −i`
  for `(3/2, 1/2), (1/2, 3/2)`, `i` otherwise ([FC99] eq 31).
- `I ⊃ C5`: phases of the multiplicity-free SAFs of lowest `l` (`l = 1, 2, 3`
  for `F1, H, F2, G`) chosen first, then the projection operators are
  required to reproduce them ([FCD99a] §3 Step 5); `I5/2` phase factors
  `ε_{μμ̄}` tabulated ([FCD99b] Table 3, eq 38c).
- `O ⊃ T ⊃ C2`: `√ϱ^{ν̄} = exp(−ν̄πi/2)` so that `G3/2` SAFs with different
  multiplicity labels agree up to an overall phase ([CFP00] eq 43b).
- `O ⊃ T ⊃ C3`: `(−1)^{ν+1/2}` so that `D(C2b)` read from the ISO agrees with
  the `j = 3/2` SAF ([F01] eq 51b).
- Dihedral groups: the universal choice `D^{E_m}(C2x) = (−1)^ξ [[0,1],[1,0]]`
  independent of `m` and `n` ("a key point in the unified solution"), with
  `E1` basis `(x, y)` for `D_n ⊃ C2x` and `(x+iy, x−iy)` for `D_n ⊃ C_n`
  ([CFMB00] eqs 2.9, 7.1).

Consequence for version 2: the eigenfunction method itself determines a
chain-adapted basis only up to one phase per basis vector (and up to the
multiplicity basis); every paper closes the gap by an external convention
(Butler agreement, real SAFs, hermiticity, associated-irrep relations). The
tool's own rule (ticket #43) has to be stated independently; the papers'
matrices are then reference data related to ours by a diagonal phase matrix
(plus the permutation of section 5.2).

---

## 2. Chains covered

Element labels: the cubic papers use Bradley & Cracknell's symbols
(`C2x, C2y, C2z, C31⁺…C34⁺, C31⁻…C34⁻, C4x^±…, C2a…C2f, I, σ_x…, S6…, S4…,
σ_da…σ_df`) with the index tables [CP99] Table 3 (`Oh`, 1–48) and Table 4
(`D6h`, 1–24) = [PC99] Tables I–II; the `T†` table numbers `E, C2x, C2y, C2z,
C31⁺..C34⁺, C31⁻..C34⁻` as 1–12 ([CF98a] Table III). The icosahedral papers
number the six `C5,j`, ten `C3,j`, fifteen `C2,j` axes by the vertices of
Fig. 1 ([FCD99a] Table 1; [FCMB99] Table I) and the 60 elements 1–60 ([FCD99a]
Table 2). Dihedral groups: `D_n = C_n ∧ C2x = {C_nz^k, C2(k) = C_nz^k C2x :
k = 0..n−1}`, the axis of `C2(k)` at azimuth `kπ/n` ([CFMB00] eq 2.1a).

Quantum numbers throughout: `C_nz = exp(−2πi J_z/n)`, eigenvalue `ρ_m =
exp(−2πim/n)`, `m` integer (s.v.) or half-integer (d.v.); `(ρ_m)^{n+k} =
θ ρ_m^k`, `C_nz^{n+k} = θ C_nz^k` ([CF98a] eqs 14–16; [FC99] eqs 5–8;
[CF98c] eq 3c lists `m` for `n = 2..5`).

| Chain | Paper, where | CSCO used | Irrep/component labels | d.v. treatment | Factor system / section |
| ----- | ------------ | --------- | ---------------------- | -------------- | ----------------------- |
| `T ⊃ D2` (s.v.) | [P89] test run 1 pp. 367–368; [CF98b] §III; [CFP00] | `C_T = Σ C3i⁺`, `C(s) = 17C + (5C2x + C2z)` (program, `D2` CSCO `5C2x + C2z`); factorization lemma `P^{(F)}_{νν'} = ℘_{νν'} P_ρ` | `A_m` (`m = 0, ±1` ↔ `A, 1E, 2E`), `F` with components `(νν') = (C2x, C2y) eigenvalues = (+−, −+, −−) ≡ (x, y, z)` | not canonical for d.v. (`C2` instead) | `T†` table [CF98a] Table III: `(C2x)² = Ẽ`, `(C31⁺)² = C̃31⁻` |
| `T ⊃ C3` | [CF98a] §IV; [CF98c] §III (CG) | `C = Σ C3i⁺`, `C(s) = C31⁺`; double cosets `T = C3(e + C2z)C3`, `q = 2` | `A_m`, `F` (weights `1, 0, −1`), `E1/2` (`±1/2`), `E3/2` (`3/2, 1/2`), `E−3/2` (`−3/2, −1/2`); eigenvalues Table I | unified; `E−3/2 = (E3/2)*` | as above; `ρ_m^3 = θ` |
| `T ⊃ C2` (d.v.) | [CF98b] §IV; [CF98c] §IV–V (CG) | `C(s) = C2z`, eigenvalue `ϱ_m = exp(−πmi) = ∓i`; `T = C2(e + C2x + C31⁻ + C31⁺)C2`, `q = 4` | `E1/2, E3/2, E−3/2` components `m = 1/2, −1/2`; `F` re-expressed in a `T ⊃ C2` real basis in [CF98c] eqs 30–31 | d.v. only (`D2` redundant) | same table |
| `O ⊃ C4` | [FC99]; [F03] (CG) | `C = Σ C4^±` (6C4), `C(s) = C4z`; `O = C4(e + C2y)C4 + C4 C2c C4`, `M_i = 4, 4, 1` | modified Mulliken `A0, A2, E (μ = 0, 2 ≡ θ, ε), F0 (0, ±1), F2 (2, ±1), E1/2 (±1/2), E3/2 (±3/2), G3/2 (±1/2, ±3/2)`; Table III | unified | `O†` table [FC99] Table II = [FCMB99] Table VI = [PC99] Table III (Herzig's table with `C2b, C2c, C2f` reversed) |
| `O ⊃ T ⊃ C3` | [F01]; [PC99] §Examples, Table XVII (CG, program) | right-induced from `T ⊃ C3` with `k̂_2 = C2b`; program `C(s) = C_T + C31⁺` | `A_σ, F_σ` (`σ = 1, 2`), `E = (A1, A−1)`, `E_ν` (`ν = 1/2, 3/2`), `G3/2 = (E3/2 ⊕ E−3/2)`; subduction [F01] Table II | unified; `G3/2` block-diagonal on `T` | same `O†` table |
| `O ⊃ T ⊃ D2 ⊃ C2` | [CFP00] | right-induced with `k̂_2 = C4z⁻`; s.v. uses `T ⊃ D2`, d.v. uses `T ⊃ C2` | `A_σ`, `E = (A1, A−1)`, `F_σ (x, y, z)`, `E_ν (±1/2)`, `G3/2` as `(E±3/2, ±1/2)`; Tables III–IV | unified | same |
| `O ⊃ D4 ⊃ D2` / `O ⊃ D4 ⊃ C2`, `O ⊃ D3 ⊃ C3`, `O ⊃ D3 ⊃ C2b`, `O ⊃ D2'` | [CFP00] §O ⊃ D4 (Tables VII–VIII, eqs 56, 59); [F01] §7A (Table IV); [PC99] Tables X–XVI | program `C(s)` of Table XVI; `O ⊃ D4 ⊃ D2 ≡ O ⊃ C4` for `E`, `≡ O ⊃ D2` for `F`; `O ⊃ D4 ⊃ C2 ≡ O ⊃ C2` for `E1/2, E3/2` | Koster labels of subgroups | — | same |
| `Td ⊃ D2d ⊃ C2`, `Td ⊃ C3v ⊃ C3`, `Td ⊃ T ⊃ …` | [CP99] §7.2 menu, test run 1 (Γ6, Γ7, Γ8 for `Td ⊃ D2d ⊃ C2z`); [CFP00] §Td (`Ŝ4z = C4z Î`, parity rule eq 55) | program | Koster | d.v. | `Td† ≅ O†` |
| `Oh`, `Th` chains | [CP99] §7.2; [PC99] eqs 14 and p. 76 (equivalences, 8 nonequivalent chains for `Oh`) | `C + 50 Î` | Koster `Γ_i^±` | trivial (`× Ci`) | inversion untouched |
| `I ⊃ C5` (s.v.) | [FCD99a] | `C = Σ_j (C5,j + C5,j⁻¹)`, `C(s) = C5z = C5,6`; `I = C5(e + C2,8)C5 + C5 C2,13 C5 + C5 C2,5 C5`, `M_i = 5, 5, 1, 1` | `A (0), F1 (0, ±1), F2 (0, ±2), G (±1, ±2), H (0, ±1, ±2)`; Table 3 | s.v. only | `I†` table "available upon request" ([FCMB99] p. 2; Euler angles Table IV) |
| `I ⊃ D5 ⊃ C5` (d.v.) | [FCD99b] | same `C`, `C(s) = C5z`, degeneracy `±5/2` lifted by `(C2,8, C̄2,8)` | `E1/2 (±1/2), E3/2 (±3/2) [= A&H E7/2, Butler 1/2~], G3/2 (±1/2, ±3/2), I5/2 (±1/2, ±3/2, [5/2±])`; Table 1 | d.v. only | same |
| `Ih ⊃ C5 × Ci` | [CP97] eq 7 (principal components only, symmetrized boson representation) | CSCO of `I` | Mulliken `A, F1, F2, G, H` with `g/u` | s.v. | — (out of scope of the map) |
| `D_n ⊃ C_n`, `D_n ⊃ C2x` (any `n`, `D∞`) | [CFMB00] §II | none (factorization lemma); irreducible-basis criterion `(J̇_z, C2x)` eq 2.22 | `A_{ρν}` / `A_{m0ν}` (`ν = ±1` s.v., `±i` d.v. for `ρ = −1`, odd `n`), `E_m`, `m = 1/2, 1, …, (n−1)/2`; Table I maps to Koster/Mulliken for `D2..D6` | unified by `ξ`, `θ` | untilded set `{C_nz^k, C_nz^k C2x : 0 ≤ k < n}`, `C_nz^{n+k} = θ C_nz^k`, `C2(2n+k) = θ C2(k)`; `(−1)^{1/2} ≡ i` |
| `C_nv ⊃ C_n`, `C_nv ⊃ C_σx` (or `C_σy`) | [CFMB00] §III | isomorphism `C2(k) → σ(k)`, `C2x → σ_x` (eq 3.2) | as `D_n` | `σ_{−n} = θ σ_n`, `σ_n = Î C2n` | as `D_n` |
| `D_nd` (even `n`) `⊃ S_2n`, `⊃ C2x` | [CFMB00] §IV | isomorphism to `D_2n`: `S'_{2n,z}^k → C_{2n,z}^{(n+1)k}`, `C2(k) → C2(2k)`, `σ(k) → C2(2k+1)` (eq 4.7) | `A_{ρν}`, `E_m`, `m = 1..n−1` or `1/2..n−1/2` | `S_nz ≡ θ σ_z C_nz` (Koster's definition, so that `Î C3z^± = S6z^∓` etc. hold for both types, eq 4.3–4.4) | as `D_2n` |
| `D_nh` (odd `n`) `⊃ C_nh`, `⊃ C2x` | [CFMB00] §V | isomorphism to `D_2n`: `C_nz^k → C_{2n,z}^{2k}`, `C2(k) → C2(2k)`, `S_nz^{(k)} = σ_z C_nz^k → C_{2n,z}^{2k+n}`, `σ(k) → C2(2k+n)` (eq 5.4); `D_nh† = C_σz† ∧ D_n†` (not a direct product in d.v.) | as `D_2n` | `S_nz^{(2k+1)} = θ^{k+1} S_nz^{2k+1}` (eq 5.3c) | as `D_2n` |
| `D_nd` (odd `n`), `D_nh` (even `n`) | [CFMB00] eq 4.1 | `D_n × Ci`, trivial | — | — | — |
| `D6 ⊃ …`, `D6h` | [CP99] Table 2; [FCMB99] Tables II, V; [PC99] Tables II, IV | `5(C6z⁺ + C6z⁻) + Σ C2''` | Koster | d.v. table `D6†` | Euler angles with `ω, θ, φ` |

Chain-adapted matrices printed for the coset generators (all exact, see the
catalogue in section 4 for pages):

- `T ⊃ C3`: `D^F(C31⁺) = diag(ρ_1, 1, ρ_{−1})`, `D^F(C2z)_{mm'} = 2/3 − δ_{mm'}`
  (every off-diagonal entry `2/3`, diagonal `−1/3`); `D^{E1/2}(C31⁺) =
  diag(ρ_{1/2}, ρ_{−1/2})`, `D^{E3/2}(C31⁺) = diag(ρ_{3/2}, ρ_{1/2})`,
  `D^{E1/2}(C2z) = D^{E3/2}(C2z) = i [[−√(1/3), √(2/3)], [√(2/3), √(1/3)]]`,
  `D^{E−3/2} = (D^{E3/2})*` ([CF98a] eqs 37–38).
- `T ⊃ D2`: `D^F(C2z) = diag(−1, −1, 1)`, `D^F(C31⁺) = [[0,0,1],[1,0,0],[0,1,0]]`
  ([CF98b] eq 17); basis relation to `T ⊃ C3`: `(ψ_x, ψ_y, ψ_z) = (1/√3)
  [[ε,1,ε*],[ε*,1,ε],[1,1,1]] (ψ_1, ψ_0, ψ_{−1})` (eq 23).
- `T ⊃ C2` (d.v.), rows `m = 1/2, −1/2`: `D(C2z) = diag(−i, i)`,
  `D^{(λ)}(C31⁺) = (λ/2√2) [[√i*, −1],[1, √i]]`, `D^{(λ)}(C2x) = [[0, i√i],
  [√i, 0]]*`, `D^{(λ)}(C31⁻) = (λ*/2√2) [[√i, 1],[−1, √i*]]`, `λ = 2, 2ε,
  2ε*` ([CF98b] eq 38); re-phased symmetric form `D(C31⁺) = (λ/2√2)
  [[√i*, −i],[−i, √i]]`, `D(C2x) = [[0, √i],[i√i, 0]]*` ([CF98c] eq 33).
- `O ⊃ C4`: `D^E(C4z) = diag(1, −1)`, `D^E(C2y) = 1`, `D^E(C2c) =
  [[−1/2, √3/2],[√3/2, 1/2]]`; `D^{F_σ}(C4z) = diag(−i, i, (−1)^{σ/2})` in
  the order `μ = 1, −1, σ`, `D^{F_σ}(C2y) = [[0,1,0],[1,0,0],[0,0,−1]]`,
  `D^{F_σ}(C2c)` with entries `±1/2`, `±1/√2`, `0` *(signs extraction
  uncertain)*; `D^{E_μ}(C4z) = diag(ρ_μ, ρ_{−μ})`, `D^{E_μ}(C2y) = [[0,−1],
  [1,0]]`, `D^{E_μ}(C2c) = −i [[1/√2, 1/√2],[1/√2, −1/√2]]` *(signs
  uncertain)*; `D^{G3/2}(C4z) = diag(ρ_{−3/2}, ρ_{−1/2}, ρ_{1/2}, ρ_{3/2})`,
  `D^{G3/2}(C2y)` an anti-diagonal `±1` matrix, `D^{G3/2}(C2c) = (i/√2)
  [[−A, B],[B, A]]` with `A, B` 2×2 blocks of entries `±1/2, ±√3/2` ([FC99]
  eqs 32–33). Reduced operator matrices for `O ⊃ C4` are restated
  algebraically in [F03] eqs 17–18 (e.g. `d^{E}_{μμ̄}(C2c) = −ρ_μ/2` on the
  diagonal, `√3/2` off it; `d^{E_ν}_{μμ}(C2c) = (−1)^{μ+ν} i/√2`,
  `d^{E_ν}_{μ,−μ}(C2c) = −i/√2`; `d^{G3/2}_{μμ}(C2c) = (−1)^{μ+|μ|+1} i/√(24)`…).
- `O ⊃ T ⊃ C3`: `D^{(A_σ)}(C2b) = (−1)^{σ+1}`, `D^{(E)}_{μμ̄}(C2b) =
  δ_{μ,−μ̄}`, `D^{(F_σ)}_{μμ̄}(C2b) = (−1)^σ δ_{μ,−μ̄}` ([F01] eq 28);
  `D^{(E_ν)}_{μμ̄}(C2b) = (−1)^{μ̄−ν} δ_{μ,−μ̄}`, `D^{(G3/2)}(C2b) =
  [[0,0,0,−1],[0,0,−1,0],[0,1,0,0],[1,0,0,0]]` in the order `(E3/2 3/2),
  (E3/2 1/2), (E−3/2 −1/2), (E−3/2 3/2)` ([F01] eq 53); the `T` elements by
  subduction (eqs 27, 52).
- `O ⊃ T ⊃ D2 ⊃ C2`: `D^{(A_σ)}(C4z⁻) = (−1)^{σ+1}`, `D^{(E)}_{νν̄}(C4z⁻)
  = δ_{ν,−ν̄}`, `D^{(F_σ)}(C4z⁻) = (−1)^σ [[0,−1,0],[1,0,0],[0,0,−1]]` in
  the order `x, y, z` ([CFP00] eqs 29–30); d.v. `D^{(E_ν)}(C4z⁻) =
  (−1)^{ν−1/2} diag(i^{1/2}, i^{−1/2})` (rows `1/2, −1/2`), `D^{(G3/2)}(C4z⁻)`
  anti-diagonal with entries `ϱ^{±1/4}`, `ϱ = exp(−πi)` ([CFP00] eqs 45–46).
- `I ⊃ C5` (s.v.): `D^{F1}(C5,6) = diag(ρ_1, ρ_0, ρ_{−1})`, `D^{F1}(β̂_2) =
  [[0,0,1],[0,−1,0],[1,0,0]]`, `D^{F2}(β̂_2) = −[[0,0,1],[0,1,0],[1,0,0]]`,
  `D^G(C5,6) = diag(ρ_2, ρ_1, ρ_{−1}, ρ_{−2})`, `D^G(β̂_2) = [[0,0,0,−1],
  [0,0,1,0],[0,1,0,0],[−1,0,0,0]]`, `D^H(C5,6) = diag(ρ_2, ρ_1, ρ_0, ρ_{−1},
  ρ_{−2})`, `D^H(β̂_2) = antidiag(1, −1, 1, −1, 1)`; `D^{F}(β̂_{3,4}) =
  (1/√5)[…]` with entries `±2c_1, ±2c_2, ±√2, ±1`, `D^G(β̂_{3,4}) =
  (1/√5)[…]` with entries `1, ±2c_k, 1/2…`, `D^H(β̂_{3,4}) = (1/5)[…]` with
  entries `4c_k, 4c_k², √6, −1` ([FCD99a] eqs 25a–d). Lowest-`l` SAFs:
  `|F1 μ⟩ = |1 μ⟩`, `|H μ⟩ = |2 μ⟩` (Table 5), i.e. `D^{F1}` and `D^H` are the
  Wigner `D^1`, `D^2` matrices in the `|l m⟩` basis with `m = μ`.
- `I ⊃ D5 ⊃ C5` (d.v.): `D^{E1/2}(C5,6) = diag(ρ_{1/2}, ρ_{−1/2})`,
  `D^{E1/2}(β̂_2) = [[0,1],[−1,0]]`, `D^{E1/2}(β̂_3) = (−i/2)
  [[s_1^{−1}, s_2^{−1}],[s_2^{−1}, −s_1^{−1}]]`, `D^{E1/2}(β̂_4) = (i/2)
  [[−s_2^{−1}, s_1^{−1}],[s_1^{−1}, s_2^{−1}]]` with `s_k = sin(2πk/5)`;
  `D^{E3/2}` analogous with `ρ_{±3/2}`; `D^{G3/2}(C5,6) = diag(ρ_{3/2},
  ρ_{1/2}, ρ_{−1/2}, ρ_{−3/2})`, `D^{G3/2}(β̂_2) = [[0,0,0,−1],[0,0,1,0],
  [0,−1,0,0],[1,0,0,0]]`, `D^{G3/2}(β̂_{3,4}) = ±(i/5)[…]` with entries
  `±t_k, ±(2/3)^{1/2}·… s_k` *(extraction uncertain)*; `D^{I5/2}(C5,6) =
  diag(ρ_{5/2}, ρ_{3/2}, ρ_{1/2}, ρ_{−1/2}, ρ_{−3/2}, ρ_{5/2})` with the two
  `−1` entries labelled `[5/2±]`, `D^{I5/2}(β̂_2)` with `∓i` in the corners
  and `±1` anti-diagonal inside, `D^{I5/2}(β̂_{3,4}) = ±(i/√5)[…]` with
  `ε_1 = 2s_2 + i s_1`, `ε_2 = (1/2)^{1/2}(−t_2 + i t_1)`, `ε_3 =
  (8/5)^{1/2}(s_2 c_2 + i s_1 c_2)` *(extraction uncertain)* ([FCD99b] eqs
  39a–e, Table 2 for `s, c, t` values).
- `D_n ⊃ C_n` (any `n`): `D^{E_m}(C_nz^k) = diag(e^{−imφ_k}, e^{imφ_k})`,
  `D^{E_m}(C2(k)) = (−1)^ξ [[0, e^{−imφ_k}],[e^{imφ_k}, 0]]`, `φ_k = 2kπ/n`
  (eq 2.12a); `D_n ⊃ C2x` real basis: `D^{E_m}(C_nz^k) = [[cos mφ_k,
  −sin mφ_k],[sin mφ_k, cos mφ_k]]`, `D^{E_m}(C2(k)) = (−1)^ξ [[cos mφ_k,
  sin mφ_k],[sin mφ_k, −cos mφ_k]]` (eq 2.12b); one-dimensional:
  `χ^{A_{ρν}}(C_nz^k) = ρ^k`, `χ(C2(k)) = ν ρ^k` (eq 2.13a); `D_nd`:
  `D^{E_m}(S_{2n,z}^k) = θ^k diag(e^{−iw_2}, e^{iw_2})`, `w_2 = (n+1)kmπ/n`,
  `D^{E_m}(σ(k)) = (−1)^ξ [[0, e^{−iw_3}],[e^{iw_3}, 0]]`, `w_3 =
  (2k+1)mπ/n` (eq 4.9); `D_nh`: `D^{E_m}(S_nz^{(k)}) = diag(e^{−iw_4},
  e^{iw_4})`, `w_4 = (2k+n)mπ/n` (eq 5.5) ([CFMB00]). Numerical 2×2 matrices
  for every element of `D4, C4v, D2d, D2, C2v` and `D6, C6v, D3h, D3, C3v`
  in the real `⊃ C2x` basis: [CFMB00] Tables VI–VIII.

The chain basis is coordinate-free: "Up to this point, the coordinate system
is irrelevant. For obtaining the SAFs, we need to specify the coordinate
system" ([CF98a] p. 471). Coordinates enter only for symmetry-adapted
functions (out of the map's scope): `T ⊃ D2` uses the Bradley–Cracknell
frame (`C2` axes along `x, y, z`, `C31⁺` toward vertex 1), `T ⊃ C3` the frame
`x'y'z'` with `z'` along `C31⁺` (`R(0, cos⁻¹(1/√3), π/4)`), `T ⊃ C2` the frame
rotated by `−π/4` about `z` to make the SAFs real, `O ⊃ T ⊃ C3` the frame
with `C2b` along `y'` and `C31⁺` along `z'` ([CF98a] Figs 1–2; [CF98b] Fig 2;
[F01] Fig 2).

Canonical chains ([P89] §4.1(e) "the subgroup has to be canonical, i.e. the
subduction is multiplicity free"; [PC99] pp. 75–76): for any point group one
intermediate subgroup suffices; the nonequivalent chains are seven for `O`,
nine for `Td`, eight for `Oh` ([PC99] p. 76 with the lists); `O ⊃ D2`
(`D2 = {E, C2x, C2y, C2z}`) is not canonical because of `E`, `O ⊃ D2'`
(`D2' = {E, C2c, C2e, C2y}`) is; `D4 ⊃ C2z` and `D2d ⊃ C2z` are canonical for
d.v. but not for s.v. ([PC99] p. 76); `D2` is redundant in `T ⊃ D2 ⊃ C2` for
d.v. and `C2` for s.v. ([CF98b] §II, §IV; [CF98c] §IV); `T` is redundant in
`O ⊃ T ⊃ C3` except for `E` and `G3/2`, `D3` in `O ⊃ D3 ⊃ C3` except for
`G3/2` ([F01] Tables II, IV); `D5` is redundant in `I ⊃ D5 ⊃ C5` except for
`I5/2` ([FCD99b] Modification 3). Equivalences `≡` (identical irreducible
bases) are listed in [PC99] p. 75: `T ⊃ D2 ⊃ C2 ≡ T ⊃ D2` (s.v.) `≡ T ⊃ C2`
(d.v.); `O ⊃ D3 ⊃ C3 ≡ O ⊃ C3`; `O ⊃ D4 ⊃ C4 ≡ O ⊃ C4`; `Oh ⊃ C4v ⊃ C4 ≡
Oh ⊃ C4`; `Th ⊃ T ⊃ C3 ≡ Ci × (T ⊃ C3)`, etc. Butler's caution that the
`D_{mn} ⊃ D_n` branch "cannot be chosen to be factorizable under all
branching schemes simultaneously" is quoted in [CFMB00] §VIII.

---

## 3. Clebsch–Gordan coefficients

**Definition and indices.** In the algebraic papers the coupled state is

```
|λ m; m̄⟩ ≡ Ψ^{(λ)m̄_1 m̄_2}_m = Ṅ P^{(λ)m̄}_m |λ_1 m̄_1, λ_2 m̄_2⟩
             = Σ_{m_1 m_2} C^{(λ)m̄, m}_{λ_1 m_1, λ_2 m_2} |λ_1 m_1⟩|λ_2 m_2⟩,
```

lower indices = product components, upper = coupled irrep, component and
multiplicity label `m̄ = (m̄_1, m̄_2)` of the trial state ([CF98c] eqs 9a–b,
10a; [F03] eq 10). Explicitly

```
C^{(λ)m̄, m}_{λ_1 m_1, λ_2 m_2} = ⟨m_1 m_2| P^{(λ)m̄}_m |m̄_1 m̄_2⟩ / ⟨m̄_1 m̄_2| P̊^{(λ)m̄}_m̄ |m̄_1 m̄_2⟩^{1/2}
                          = sqrt(h_λ/|G|) Ṅ Σ_a D^{(λ_1)}_{m_1 m̄_1}(R_a) D^{(λ_2)}_{m_2 m̄_2}(R_a) D^{(λ)}_{m m̄}(R_a)*
```

([CF98c] eq 9b, attributed to Dirl 1979), and with the reduced operator

```
C^{(λ)m̄, m}_{λ_1 m_1, λ_2 m_2} = N Σ_{i=1}^{q} M_i^{−1} d^{(λ_1)}_{m_1 m̄_1}(β_i) d^{(λ_2)}_{m_2 m̄_2}(β_i) d^{(λ)}_{m m̄}(β_i)*,
m ≐ m_1 + m_2,   m̄ ≐ m̄_1 + m̄_2                                              (CF98c eq 11; F03 eq 11)
```

the "magnetic quantum number conservation rule" (congruence mod `n`). The
programs use `ψ^{(λ'')τ}_{a''} = Σ C^{λ a, λ' a'}_{λ'' a'' τ} ψ^λ_a ψ^{λ'}_{a'}`
and print rows `(Γ, a, MUL)` against columns `a a'` ([CP99] eq 34, §8.3;
[P89] §5 table heading).

**Normalisation and phase.** The coupled state is normalised and `N`, `Ṅ`
are real positive ([CF98c] eqs 9a, 10b–c), so the coefficient of the trial
product state `|λ_1 m̄_1, λ_2 m̄_2⟩` in the component `m = m̄` is **real and
positive**; all other components of the same coupled irrep follow from the
ISO (equivalently from the irrep matrices), exactly as in the programs' rule
"the same relative phase as for the irreps" ([CP99] p. 82). The only
convention is therefore the choice of trial state per `(λ_1, λ_2, λ)`:
[CF98c] Table II lists them for `T ⊃ C3` (e.g. `F × E1/2 → E1/2: |0, 1/2⟩`;
`E1/2 × E1/2 → A0: |1/2, −1/2⟩`, `→ F: |1/2, −1/2⟩`; `E3/2 × E3/2 → A−1:
|1/2, 3/2⟩`, `→ F: |3/2, 3/2⟩`); for `T ⊃ D2 ⊃ C2` the trial states are
`|A_m, F_z⟩`, `|F_x, F_x⟩`, `|F_x, F_y⟩ ± |F_y, F_x⟩`, `|A_m, E_ν 1/2⟩`,
`|F_z, λ 1/2⟩`, `|λ_1 1/2, λ_2 −1/2⟩` ([CF98c] §IV, eqs 27–29 with the
selection rules `ν_1 ν_2 = ν`, `ν_1' ν_2' = ν'` on the `D2` quantum numbers
and `ς̄_1 + ς̄_2 = 1/2` on the `C2` ones); for `O ⊃ C4` they are implicit in
[F03] eqs 20–22 (e.g. `F_{σ_1} × F_{σ_2} → A_μ`: `(1/√3)[|σ_1, σ_2⟩ −
|1, μ−1⟩ − |−1, μ+1⟩]`, `E_{ν_1} × E_{ν_2} → A_μ`: `(1/√2)[|μ̄_1, μ̄_2⟩ −
ρ_μ |−μ̄_1, −μ̄_2⟩]`). For `D_n` ([CFMB00] eqs 2.24–2.26) the trial state is
`|m_1, m_2⟩` with `m_1 > 0`:

```
A_{m0 ν} × E_{m2} → E:  C^E_m = δ̇_{m, m0+m2} |A, m2⟩ + ν (−1)^{ξ2−ξ} δ̇_{m, −(m0+m2)} |A, −m2⟩,
E_{m1} × E_{|m2|} → A_{m0 ν}: C = (1/√2)(|m1 m2⟩ + ν* (−1)^{ξ1+ξ2} |−m1, −m2⟩),  m1 + m2 ≐ m0 ∈ {0, n/2},
E_{m1} × E_{|m2|} → E_m:      C^E_m = δ̇_{m, m1+m2} |m1 m2⟩ + δ̇_{m, −(m1+m2)} (−1)^{ξ1+ξ2−ξ} |−m1, −m2⟩,
```

with `δ̇` the Kronecker delta mod `n`. The `D3` example is written out in
eqs 2.27a–d.

**Coupling multiplicity.** For the point groups the s.v. products are
multiplicity-free except `F × F → 2F` in `T` ([CF98c] eq 19a), `F_σ × G3/2 →
2G3/2` and `G3/2 × G3/2 → 2(F0 + F2)` in `O` ([F03] eq 19), `Γ8 × Γ8 → 2Γ4 +
2Γ5` in [PC99] Table XVII; `(ijk) ≤ 1` for s.v. reps, `> 1` possible for d.v.
([PC99] p. 77). The algebraic papers resolve it by the **intrinsic quantum
number** `m̄ = (m̄_1, m̄_2)`: two different trial states give two coupled
copies which "are in general not orthogonal with respect to `m̄`" ([CF98c]
p. 5522, 5525: `F × F → F` from `|F0, F0⟩` and `|F1, F−1⟩`, orthogonalised
before printing Table III; [F03] Table II footnote: `Ψ(G3/2)^1 ≡
Ψ^{(G3/2) 1, −1/2}`, `Ψ(G3/2)^2 ≡ Ψ^{(G3/2) −1, 3/2}`; Table III footnote:
`Ψ(F0)^1 ≡ Ψ^{1/2, −1/2}`, `Ψ(F0)^2 ≡ Ψ^{3/2, −3/2}`, `Ψ(F2)^1 ≡
Ψ^{1/2, 3/2}`, `Ψ(F2)^2 ≡ Ψ^{3/2, 1/2}`). Orthogonalising introduces new
square roots: the second `G3/2` copy of `F0 × G3/2` in [F03] Table II has
entries `√(5/6), 1/√15, 1/√10, √(3/10), √(3/5)` while every multiplicity-free
`O ⊃ C4` coefficient is in `{±1, ±1/2, ±√3/2, ±1/√2, ±1/√3, ±√(2/3), ±1/√6}`.
No symmetric/antisymmetric rule is used for `λ_1 = λ_2`; the symmetry of the
copies is reported afterwards: for `D_n`, "`E_{2m}` and `A_{n/2±}` are in the
symmetrized part of `E_m × E_m`, while `A_{01}` (`A_{02}`) is in the symmetric
part for integer (half-integer) `m`" ([CFMB00] after eq 2.26c); for `T ⊃ D2`,
`C^{(F)}_{z, r} = (1/√2)[|F_x F_y⟩ + (−1)^r |F_y F_x⟩]`, `r = 0, 1` ([CF98c]
eq 28b) — here the two copies happen to be the symmetric and antisymmetric
ones. The programs leave the multiplicity basis to the eigensolver (section
1.3).

**CG series.** Obtained from the weights (additive `m` mod `n`) by the
subtraction method of Lie groups rather than from characters ([CF98c] §III.A:
`A_{m1} × A_{m2} = A_{m1+m2}`, `A_m × F = F`, `F × F = A0 + A1 + A−1 + 2F`,
`F × E_ν = E1/2 + E3/2 + E−3/2`, `A_{±1} × (E1/2, E±3/2, E∓3/2) = (E±3/2,
E∓3/2, E1/2)`, `E_{ν1} × E_{ν2} = A_m + F` with `ρ_m = λ_1 λ_2 / 4`; [F03]
eq 19 for `O`, e.g. `A_{μ1} × A_{μ2} = A_μ`, `μ = μ_1 + ρ_{μ1} μ_2`, `E ×
E = A0 + A2 + E`, `E × F_σ = F0 + F2`, `F_{σ1} × F_{σ2} = A_μ + E + F0 + F2`,
`E × E_ν = G3/2`, `E × G3/2 = E1/2 + E3/2 + G3/2`, `F_σ × E_ν = E_{ν'} +
G3/2`, `E_{ν1} × E_{ν2} = A_μ + F_σ` with `μ = σ = ν_1 + (−1)^{ν1+ν2} ν_2`,
`E_ν × G3/2 = E + F0 + F2`, `G3/2 × G3/2 = A0 + A2 + E + 2(F0 + F2)`; [CFMB00]
eq 2.23 for `D_n`: `E_{m1} × E_{m2} = A01 + A02 + E_{2m1}` (`m_1 = m_2`), `=
A01 + A02 + A_{n/2,1} + A_{n/2,2}` (`m_1 = m_2 = n/4`), `= A_{n/2,1} +
A_{n/2,2} + E_{m1−m2}` (`m_1 + m_2 = n/2`), `= E_{m1+m2} + E_{|m1−m2|}`
otherwise).

**Symmetry relations.** Exchange of the factors: `C^{Γ_2 m_2, Γ_1 m_1}_{Γ r m}
= (−1)^{Γ_1+Γ_2+Γ+r} C^{Γ_1 m_1, Γ_2 m_2}_{Γ r m}` with `(−1)^Γ = 1, −1, i`
for `A_m, F, E_ν` of `T` and `r = 0, 1` the multiplicity label of `F × F → F`
([CF98c] eq 23, quoting Butler & Wybourne). Sign reversal of all components
through a binary rotation: for `O ⊃ C4`, `C^{(λ)μ̄_1μ̄_2, −μ}_{λ_1 −μ_1, λ_2 −μ_2}
= D^{(λ)−1}_{−μ μ}(C2y) D^{(λ_1)}_{−μ_1 μ_1}(C2y) D^{(λ_2)}_{−μ_2 μ_2}(C2y)
C^{(λ)μ̄_1μ̄_2, μ}_{λ_1 μ_1, λ_2 μ_2}` ([F03] eq 14). Reality: all `T ⊃ C3`,
`T ⊃ C2` and `O ⊃ C4` coefficients are real ([CF98c] p. 5525, 5533; [F03]
p. 549); `T ⊃ D2 ⊃ C2` coefficients are complex "in agreement with Butler
that some `3jm` factors of `T ⊃ D2` must be complex" ([CF98c] p. 5530);
`O ⊃ T ⊃ C3` coefficients involve `exp(±iπ/6)`, `exp(iπ/3)` ([PC99] Table
XVII). Isomorphic groups: "with suitable choices of phase the CG coefficients
of an improper dihedral group are identical to those of the dihedral group it
is isomorphic to" (Butler), so [CFMB00] gives CG coefficients for `D_n` only
and maps `C_nv → D_n`, `D_nd, D_nh → D_2n` (§II.G end, §VIII).

---

## 4. Printed tables

"Exact" means closed forms in surds and roots of unity (the programs print
`sqrt(p/q) exp(iπ n/m)` or `sqrt(p/q) exp(i cos⁻¹(m/n))`, [CP99] abstract;
[PC99] abstract). SAF tables (symmetry-adapted functions, out of the map's
scope) are listed for completeness because they are the only place where some
phase choices are visible. Pages are those printed on the PDF.

| Paper | Table / eq | Page | Chain | Contents | Exact / decimal | Double-valued? |
| ----- | ---------- | ---- | ----- | -------- | --------------- | -------------- |
| [P89] | Test run 1 | 367–368 | `T ⊃ D2` | class operators, CSCO-I/II/III of `T`, character table with eigenvalues, s.v. irrep matrices of all 12 elements (`A, 1E, 2E, F`) | exact (`E = −exp(−iπ/3)`) | no |
| [P89] | Test run 2–3 | 368–373 | space groups `T⁴`, `Oh¹⁰` | little co-group irreps and CG coefficients of space groups | exact | s.v. shown |
| [CP97] | eq 7, Tables 1–3, cases 1–6 | 405–413 | `Ih ⊃ C5 × Ci` | principal components of `Ih` irreducible bases in the symmetrized boson representation (B12H12, C60) | exact | no (out of scope) |
| [CF98a] | Table I | 469 | `T ⊃ C3` | eigenvalues `λ`, weights `m`, `ρ_m` ↔ Koster `Γ1–Γ7` and modified Mulliken | exact | yes |
| [CF98a] | Table II | 473 | `T†` | Euler angles of the 12 elements (Fig. 1 frame) | exact | yes |
| [CF98a] | Table III | 473 | `T†` | 12×12 multiplication table with tildes | — | yes |
| [CF98a] | eqs 37–38, `(j,k)` list | 479 | `T ⊃ C3` | irrep matrices of `C31⁺`, `C2z` for `F, E1/2, E3/2` (`E−3/2` conjugate); `(j,k)` of every element in `C31^j C2z C31^k` | exact | yes |
| [CF98a] | Table IV | 483–484 | `T ⊃ C3` | s.v. SAFs `l = 1–3`, part of `6, 8` | exact | no |
| [CF98a] | Table V | 486–487 | `T ⊃ C3` | d.v. SAFs `j = 1/2–9/2` | exact | yes |
| [CF98b] | Table I | 5506 | `T†` | multiplication table (repeat) | — | yes |
| [CF98b] | eq 17, eq 23 | 5507, 5510 | `T ⊃ D2` | `D^F(C2z)`, `D^F(C31⁺)`; transformation to the `T ⊃ C3` basis | exact | no |
| [CF98b] | Table II | 5510 | `T ⊃ D2` | s.v. SAFs `l = 1–3`, part of `6, 8, 13` (B–C real harmonics) | exact | no |
| [CF98b] | eq 38, element list | 5514 | `T ⊃ C2` | d.v. matrices of `C2z, C31⁺, C2x, C31⁻` for `E1/2, E3/2, E−3/2`; `C2z^j C2x`, `C2z^j C31^± C2z^k` lists | exact (`√i`) | yes |
| [CF98b] | Table III | 5516–5517 | `T ⊃ C2` | d.v. SAFs `j = 1/2–9/2` | exact | yes |
| [CF98c] | Table I | 5524 | `T ⊃ C3` | eigenvalues and weights (repeat) | exact | yes |
| [CF98c] | eq 17 | 5524 | `T ⊃ C3` | `d(C2z)` for all irreps (algebraic) | exact | yes |
| [CF98c] | Table II | 5526 | `T ⊃ C3` | trial states for the CG coefficients | — | yes |
| [CF98c] | Table III | 5527 | `T ⊃ C3` | CG coefficients, all products `Γ1 ≤ Γ2` (s.v.×s.v., s.v.×d.v., d.v.×d.v.), `F×F→F` two orthogonalised copies | exact, real | yes |
| [CF98c] | eqs 24d, 26, 30–31, 33 | 5528–5532 | `T ⊃ D2`, `T ⊃ C2` | `D^F(C31⁺)`, d.v. `D(C31^±), D(C2x), D(C2y)`; `F` in the `T ⊃ C2` real basis | exact | yes |
| [CF98c] | eqs 21a, 28a–b, 29a–d; Table IV | 5525, 5529–5531 | `T ⊃ D2 ⊃ C2` | CG coefficients (`A×A`, `A×F`, `F×F` in eqs; `F×E_ν`, `E×E` tabulated) | exact, complex | yes |
| [CF98c] | eqs 34–35; Table V | 5532–5533 | `T ⊃ C2` | CG coefficients `F×F→A`, `F×E_ν`, `E×E` (real basis) | exact, real | yes |
| [CP99] | Tables 1–2 | 75 | `O†`, `D6†` | Euler angles (Rose), B–C labels | exact | yes |
| [CP99] | Tables 3–4 | 76 | `Oh`, `D6h` | element index 1–48, 1–24 | — | — |
| [CP99] | Tables 5–6 | 84 | subgroups and chains of `T, Td, O, Oh` | elements and CSCO of `O`, `D6` and subgroups; `C(s)` of 24 preset chains with element indices | exact | both |
| [CP99] | §7.2 menus | 82–83 | `T, Td, O, Oh` | preset chain lists (e.g. `O`: `D4⊃C2x, D3⊃C2b, D3⊃C3, T⊃C3, C4z, D4⊃D2, T⊃D2, D2'`, `T⊃C2x` d.v.) | — | both |
| [CP99] | Test run 1 | 89–90 | `Td ⊃ D2d ⊃ C2z` | character table and d.v. irrep matrices `Γ6, Γ7, Γ8` on all 24 elements, component labels `[Γ6 γ3, Γ6 γ4, …]` | exact (`E = √(1/2) e^{−iπ/4}`, `G = √(3/8)`, …) | yes |
| [CP99] | Test run 2 | 91 | `O ⊃ D4 ⊃ C2x` (s.v.) × `O ⊃ C4` (d.v.) | CG coefficients `Γ2×Γ6→Γ7`, `Γ5×Γ8→Γ6+Γ7+2Γ8` (partial) | exact | yes |
| [CP99] | Test runs 3–4 | 92–93 | space group `Oh⁷` | `D4h ⊃ C4` d.v. little co-group irreps, CG | exact | yes |
| [FC99] | Table I | 8 | `O†` | Euler angles | exact | yes |
| [FC99] | Table II | 9 | `O†` | 24×24 multiplication table with tildes | — | yes |
| [FC99] | Table III | 10 | `O ⊃ C4` | eigenvalues `λ`, weights `μ`, modified Mulliken / Mulliken / Koster | exact | yes |
| [FC99] | eqs 30–31 | 13 | `O ⊃ C4` | reduced projection operators (algebraic, with phase factors `θ`) | exact | yes |
| [FC99] | eqs 32–33, element list | 17–18 | `O ⊃ C4` | irrep matrices of `C4z, C2y, C2c` for `E, F0, F2, E1/2, E3/2, G3/2`; `(j,k)` of every element | exact | yes |
| [FC99] | Table IV | 14–15 | `O ⊃ C4` | s.v. SAFs `l = 1, 2, 3, 5, 6, 8, 12, 15` | exact | no |
| [FC99] | Table V | 16–17 | `O ⊃ C4` | d.v. SAFs `j = 1/2–9/2, 13/2, 17/2` | exact | yes |
| [FCMB99] | Table I | 5 | `I` | rotation axes by vertex labels | — | — |
| [FCMB99] | Tables II–IV | 5–7 | `D6†`, `O†`, `I†` | unique Euler angles with `ω`, `(θ, φ)` of every element (12, 24, 60) | exact | yes |
| [FCMB99] | Tables V–VI | 8 | `D6†`, `O†` | multiplication tables with tildes | — | yes |
| [FCMB99] | eq 18 | 6 | `D_n ⊃ C_n` | `E1/2` matrices of `C_nz` and `C2^{(k)}` (any `n`) | exact | yes |
| [PC99] | Tables I–II | 69 | `Oh`, `D6h` | element index | — | — |
| [PC99] | Tables III–IV | 70–71 | `O†`, `D6†` | multiplication tables | — | yes |
| [PC99] | Table V | 71 | subgroups of `O`, `D6` | elements and CSCO (`C31, D2, D2', C2v, C2v', D3, C3v, D4, C4v, D2d, T, O, D6`) | exact | both |
| [PC99] | Tables VI–XV | 72–74 | `D2, T⊃D2, T⊃C2z, T⊃C31, O⊃D2', O⊃C4z, D4⊃D2, D4⊃C2x, D2d⊃C2v', D2d⊃Cσa, D3⊃C2b, D3⊃C31` | eigenvalues of CSCO and `C(s)` ↔ Koster labels `(Γ_i, γ_j)` | exact | both |
| [PC99] | Table XVI | 76 | chains of `Td, O, Oh` | `C(s)` of 19 chains with element indices | exact | both |
| [PC99] | §Examples | 77 | `O ⊃ T ⊃ C31` | d.v. matrices `Γ6, Γ7, Γ8` of `C2x, C31⁺, C4x⁺, C2a` (`a = −1/√3`, `b = e^{−iπ/3}/√3`, `c = e^{iπ/6}/√3`, `d = √(2/3)`, `f = −√(2/3) e^{−iπ/3}` *(extraction uncertain)*) | exact | yes |
| [PC99] | Table XVII | 78 | `O ⊃ T ⊃ C31` | CG coefficients `Γ6×Γ6, Γ6×Γ7, Γ7×Γ7, Γ6×Γ8, Γ7×Γ8, Γ8×Γ8` (d.v.×d.v.→s.v.), with multiplicity | exact (`e^{±iπ/6}`, `e^{iπ/3}`) | yes |
| [FCD99a] | Tables 1–2 | 872–873 | `I†` | rotation axes; Euler angles of 60 elements | exact | yes |
| [FCD99a] | Table 3 | 874 | `I ⊃ C5` | eigenvalues `(λ, μ)` of `(C, C5)` for `A, F1, F2, G, H` | exact | no |
| [FCD99a] | eqs 24a–e | 876 | `I ⊃ C5` | reduced projection operators (algebraic) | exact | no |
| [FCD99a] | eqs 25a–d, Table 4 | 877–878 | `I ⊃ C5` | irrep matrices of `C5,6, β̂_2, β̂_3, β̂_4` for `F1, F2, G, H`; `(j,k)` of all 60 elements | exact (`c_k, √5, √2, √6`) | no |
| [FCD99a] | Table 5 | 880 | `I ⊃ C5` | s.v. SAFs `l = 1, 2, 3, 6, 8, 9, 11` | exact | no |
| [FCD99a] | eqs 37, 39 | 882 | `I ⊃ C5` | phase relations to [CP97] bases; `U` matrices to McLellan's `SO3 ⊃ I` basis | exact | no |
| [FCD99b] | Table 1 | 1050 | `I ⊃ D5 ⊃ C5` | subduction and eigenvalues `(λ, μ)` for `E1/2, E3/2, G3/2, I5/2`, with A&H and Butler names | exact | yes |
| [FCD99b] | Table 2 | 1053 | — | values of `s_k, c_k, t_k` | exact | — |
| [FCD99b] | Table 3 | 1055 | `I ⊃ D5 ⊃ C5` | phase factors `ε_{μμ̄}` of `I5/2` | exact | yes |
| [FCD99b] | eqs 38a–c | 1054–1055 | `I ⊃ D5 ⊃ C5` | reduced projection operators | exact | yes |
| [FCD99b] | eqs 39a–e, Table 4 | 1055–1056 | `I ⊃ D5 ⊃ C5` | irrep matrices of `C5,6, β̂_2, β̂_3, β̂_4` for `E1/2, E3/2, G3/2, I5/2`; `(j,k)` of all elements | exact | yes |
| [FCD99b] | Table 5 | 1057 | `I` | multiplicities of irreps of `I` in `D^j ↓ I`, `j ≤ 14`, `≤ 29/2` | integers | both |
| [FCD99b] | Table 6 | 1058 | `I ⊃ D5 ⊃ C5` | d.v. SAFs `j = 1/2–7/2, 15/2, 17/2` | exact | yes |
| [FCD99b] | eq 47 | 1059 | `I ⊃ D5 ⊃ C5` | `U` matrices to McLellan's basis (`U^{E1/2} = diag(1, i)`, …) | exact | yes |
| [CFP00] | Tables I–II | 587 | `O`, `T` | eigenvalues ↔ modified Mulliken / Mulliken / Koster | exact | yes |
| [CFP00] | Tables III–IV | 588–589 | `O ⊃ T ⊃ D2`, `O ⊃ T ⊃ C2`, `O ⊃ D4 ⊃ C2` | subduction rules with component indices | — | both |
| [CFP00] | eqs 28–30 | 591 | `O ⊃ T ⊃ D2` | s.v. matrices of `C2z, C31⁺` (from `T`) and `C4z⁻` for `A_σ, E, F_σ` | exact | no |
| [CFP00] | Table V | 593 | `O ⊃ T ⊃ D2` | s.v. SAFs `l = 1–3`, part of `6, 8, 13` | exact | no |
| [CFP00] | eqs 44–46 | 594 | `O ⊃ T ⊃ C2` | d.v. matrices of `C2z, C31^±, C2x` (from `T`) and `C4z⁻` for `E_ν, G3/2` | exact | yes |
| [CFP00] | Table VI | 596 | `O ⊃ T ⊃ C2` | d.v. SAFs `j = 1/2–9/2` | exact | yes |
| [CFP00] | Tables VII–VIII, eqs 56, 59 | 597–598 | `O ⊃ D4 ⊃ C4`, `O ⊃ D4 ⊃ D2`, `O ⊃ D4 ⊃ C2` | subduction; `C4` ↔ `(D4, C2)` labels; basis relations `E` and `G3/2` between `O ⊃ T` and `O ⊃ D4` chains | exact | yes |
| [CFP00] | eq 55 | 597 | `Td ⊃ T ⊃ D2 ⊃ C2` | `Td` from `O`: identical for even parity, label swaps `A1↔A2, F1↔F2, E1/2↔E3/2`, sign changes for `E, G3/2` for odd parity | — | yes |
| [CFMB00] | Table I | 8199 | `D2–D6` | Koster / Mulliken / modified Mulliken irrep names | — | yes |
| [CFMB00] | eqs 2.12a–b, 2.13 | 8201 | `D_n ⊃ C_n`, `D_n ⊃ C2x` | irrep matrices of every element and characters, any `n`, s.v. and d.v. | exact (algebraic) | yes |
| [CFMB00] | eqs 2.17, 2.19, 2.20 | 8202–8203 | `D_n ⊃ C_n`, `D_n ⊃ C2x`; `D3`, `C3v` | SAFs (algebraic) and lowest-`l` examples | exact | yes |
| [CFMB00] | eqs 2.23, 2.26, 2.27 | 8204–8205 | `D_n ⊃ C_n` | CG series and CG coefficients (algebraic), `D3` example | exact | yes |
| [CFMB00] | eqs 2.30–2.33 | 8206–8207 | `D∞ ⊃ C∞` | matrices, characters, SAFs, CG | exact | yes |
| [CFMB00] | eqs 3.5–3.9 | 8208–8210 | `C_nv ⊃ C_n`, `⊃ C_σx` / `C_σy`; `C∞v` | SAFs with parity dependence | exact | yes |
| [CFMB00] | eqs 4.9–4.13 | 8212–8213 | `D_nd ⊃ S_2n`, `⊃ C2x` (even `n`) | matrices, characters, SAFs | exact | yes |
| [CFMB00] | Tables II–IV | 8214 | `D2d`, `D3h` | one-dimensional SAF `m` ranges; `D2d ⊃ S4` two-dim SAF first components; minimum set of one-dim SAFs | exact | yes |
| [CFMB00] | eqs 5.5–5.6, Table V | 8215–8216 | `D_nh ⊃ C_nh`, `⊃ C2x` (odd `n`); `D3h ⊃ C3h` | matrices, characters; SAF first components | exact | yes |
| [CFMB00] | eqs 6.5–6.6 | 8217 | `O3 ⊃ D_n` | subduction frequencies (algebraic) | exact | yes |
| [CFMB00] | Tables VI–VIII | 8218–8219 | `D4, C4v, D2d, D2, C2v ⊃ C2x`; `D6, C6v, D3h, D3, C3v ⊃ C2x` | numerical 2×2 matrices of all elements, s.v. and d.v., real basis; the constituent matrices `e, r, l, k, a, b, m, n, â, b̂, m̂, n̂` | exact | yes |
| [F01] | Tables I–III | 262 | `T`, `O ⊃ T ⊃ C3`, `O` | eigenvalues ↔ labels; subduction | exact | yes |
| [F01] | Table IV | 265 | `O ⊃ D3 ⊃ C3` | subduction (`G3/2 → E1/2 + A3/2 + A−3/2`) | — | yes |
| [F01] | eqs 26–28 | 263 | `O ⊃ T ⊃ C3` | s.v. projection operators and `D(C2b)` for `A_σ, E, F_σ` | exact | no |
| [F01] | eqs 51–53 | 266 | `O ⊃ T ⊃ C3` | d.v. projection operators and `D(C2b)` for `E_ν, G3/2` | exact | yes |
| [F01] | Table V | 268 | `O ⊃ T ⊃ C3` | s.v. SAFs `l = 1–3, 6, 8` (as symmetrised `T ⊃ C3` SAFs) | exact | no |
| [F01] | Table VI | 269 | `O ⊃ T ⊃ C3` | d.v. SAFs `j = 1/2–9/2` | exact | yes |
| [F01] | §7A | 267–268 | `O ⊃ D3 ⊃ C3` | relation of `G3/2` bases (`A±3/2` of `D3` ↔ `E±3/2` of `T`, eigenvalue `∓i` of `C2b`) | exact | yes |
| [F03] | eqs 15–18 | 543–544 | `O ⊃ C4` | reduced projection operators and `d(C2y), d(C2c)` for all irreps (algebraic) | exact | yes |
| [F03] | eq 19 | 544 | `O` | CG series in modified Mulliken notation | — | yes |
| [F03] | eqs 20–22 | 544–550 | `O ⊃ C4` | CG coefficients (algebraic) for every product class | exact | yes |
| [F03] | Table I | 545 | `O ⊃ C4` | CG coefficients s.v.×s.v.→s.v. (`E×E, E×F0, E×F2, F0×F0, F0×F2, F2×F2`) | exact, real | no |
| [F03] | Table II | 546–547 | `O ⊃ C4` | CG coefficients s.v.×d.v.→d.v. (`E×G3/2, F0×E1/2, F0×E3/2, F0×G3/2, F2×E1/2, F2×E3/2, F2×G3/2`) incl. two `G3/2` copies | exact, real | yes |
| [F03] | Table III | 548–549 | `O ⊃ C4` | CG coefficients d.v.×d.v.→s.v. (`E1/2×E1/2, E1/2×E3/2, E1/2×G3/2, E3/2×E3/2, E3/2×G3/2, G3/2×G3/2`) incl. two `F0`, `F2` copies | exact, real | yes |

Not printed anywhere: complete irrep matrices for all elements of `O`
(except through `(j,k)` recipes), the `I†` multiplication table ("available
upon request", [FCMB99] p. 2, [FCD99a] p. 873), CG coefficients of `I`, CG
coefficients of `O ⊃ T ⊃ D2 ⊃ C2` and `O ⊃ T ⊃ C3` beyond [PC99] Table XVII,
any table for `n > 6` dihedral groups (covered by the formulas), and any
`Th`, `Oh`, `Ih` table (trivial by `× Ci`).

---

## 5. Differences from this tool (version 1 conventions)

### 5.1 Section and factor system

- **Model.** The rep group `{R_i, −R_i}` with `R̃_i = −R_i` is our `G†` read
  through a projective representation; `η(i, j) = ±1` of [P89] eq 2.6 is our
  `ω`. The inversion is treated as commuting and trivially lifted (`G × Ci`,
  `C + 50Î`; `Î C_n^± = S^∓` for both types), exactly as in
  `docs/spec/double-groups.md` ("the inversion lifts trivially"). Koster's
  `S_nz ≡ θ σ_z C_nz` ([CFMB00] eqs 4.2–4.5) is the consequence: the untilded
  `S_n` is `Î` times the lifted rotation part, which is our `D(i·r) = ±D(r)`
  rule. The naive `σ_z C_nz` would be the barred element in d.v.
- **Rotations by angle.** [FCMB99] takes the untilded element as the rotation
  by `ω ∈ [0, 4π)` with the Euler angles of Tables II–IV; for the printed
  `C_n^-` elements `ω = 10π/3` (`C3^-`), `11π/3` (`C6^-`), i.e. `−2π/n mod
  4π`, so the untilded `C_n^-` is the inverse of the untilded `C_n^+`
  (`C_n^{−k} ≡ (C_n^k)^{−1}`, e.g. `C3^- ≡ C̃3^2`, [FCMB99] p. 5) — the same
  as our canonical
  section (angle in `(−π, π]`). The `T†` table confirms `(C31⁺)² = C̃31⁻` and
  `(C2x)² = Ẽ` ([CF98a] Table III), our `ω(C3, C3) = −1`, `ω(C2, C2) = −1`.
- **Powers in the dihedral paper.** [CFMB00] takes the untilded set
  `{C_nz^k : 0 ≤ k < n}` with `D^{E1/2}(C_nz^k) = diag(e^{−ikπ/n}, e^{ikπ/n})`
  (eq 2.12a), i.e. the rotation angle in `[0, 2π)`, and `C2(k) = C_nz^k C2x`
  with azimuth `kπ/n ∈ [0, π)`. For `n/2 < k < n` these are the **barred**
  elements of our section (our residue rule reduces `k` to `(−n/2, n/2]`,
  and the pole of our lifted `C2'(k)` has azimuth in `(−π/2, π/2]`; e.g. for
  `D3` their `C2(2)` is the quaternion `(0, cos 2π/3, sin 2π/3, 0)` whose
  first non-zero coordinate is negative). The relative labelling is
  `ε(C_n^k) = ε(C2(k)) = −1` for `n/2 < k < n`, `+1` otherwise, and the
  papers' factor system is `ω_ref = ε ε ε ω` as in the "Section overrides"
  subsection of `docs/spec/double-groups.md`.
- **Binary axes of `O`.** The `O†` table used by [CP99], [FC99], [FCMB99],
  [PC99], [CFP00], [F01], [F03] is Herzig's (Altmann's quaternion
  parameterisation) "except that the three two-fold axes `C2b, C2c, C2f` are
  reversed" ([CP99] p. 75; [FCMB99] Remark 2; [PC99] p. 68). Their poles
  ([FCMB99] Table III, `(θ, φ)` column, *extraction partly uncertain*):
  `C2a (π/2, π/4)`, `C2b (π/2, 3π/4)`, `C2c (π/4, 0)`, `C2d (π/4, π/2)`,
  `C2e (π/4, π)`, `C2f (3π/4 or π/4, π/2 or 3π/2)`. `C2b` (`(−1, 1, 0)/√2`)
  and `C2e` (`(−1, 0, 1)/√2`) lie in our negative hemisphere, `C2a, C2c,
  C2d` in the positive one. Version 1 records the Altmann & Herzig swaps for
  `O` as `C2'(1), C2'(3), C2'(5)` in canonical labels; the B–C ↔ canonical
  label map and the resulting `ε` for the papers' table must be settled by
  the verification ticket (#48) from the multiplication table itself
  ([FC99] Table II), not from the pole table.
- **Icosahedral section.** Only Euler angles are printed ([FCMB99] Table IV;
  [FCD99a] Table 2); `C5,j^{−1}` elements carry `ω = 18π/5 ≡ −2π/5` and the
  `C2,j` elements `ω = π` with poles toward the midpoints of edges listed in
  Table 1 — the same remark as for `O` applies.
- **`D_nh` with odd `n`.** `D_nh† = C_σz† ∧ D_n†` is not a direct product
  ([CFMB00] eq 5.2); our model gets the same group through the type III
  construction. The isomorphism `D_nh ≅ D_2n` of eq 5.4 (`S_nz^{(k)} →
  C_{2n,z}^{2k+n}`) matches our `N = 2n` principal element `S_n^{−1}` of order
  `2n` for `D_nh`/`C_nh` with odd `n` (`docs/spec/irrep-labels.md`, named
  elements).

### 5.2 Quantum numbers, component order and the canonical basis

- **Sign of the quantum number.** Chen's `m` is the eigenvalue label of the
  active rotation `C_nz = exp(−2πi J_z/n)` on `|m⟩`: `ρ_m = e^{−2πim/n}`. Our
  canonical basis lists eigenvalues `ζ_N^e` with the exponent `e` descending
  and the component-`1` member of a complex pair carrying `ζ_N^{+k}`. Hence
  **`m = −e`**: Chen's component `m = +1/2` is our `e = −1/2`, his `A_1` of
  `T` (`χ(C31⁺) = ε = ζ_3^{−1}`) is our `2E`, his `A_{−1}` our `1E`, his
  `E3/2` of `T` (`χ(C31⁺) = ε`) our `2E3/2` and his `E−3/2` our `1E3/2`; for
  `D_n` with odd `n` his `A_{n/2, i}` (`D(C2x) = i`, `≡ A_{n/2,2}`) is our
  `1E_{N/2}` (`χ(s(C2'(1))) = +i`) and `A_{n/2,−i}` our `2E_{N/2}`. The
  active convention is stated explicitly: "our character table is the complex
  conjugate of those using the passive convention [Koster]" ([CP99] §9).
- **Component order.** With `m` descending from `+|m|` ([CF98a] eq 38;
  [CFMB00] eq 2.12a; [FCD99b] eqs 39) the papers' order is the **reverse** of
  our descending-exponent order for every irrep whose components are labelled
  by the cyclic quantum number alone (`E_m`, `E_{k/2}`, `F`, `G`, `H`,
  `G3/2`, `E1/2`, `E3/2`). Exceptions: `O ⊃ C4` lists `F_σ` as `(1, −1, σ)`
  ([FC99] eq 32) and `G3/2` as `(−3/2, −1/2, 1/2, 3/2)` ([FC99] eq 33b, our
  order); `G3/2` of `O ⊃ T ⊃ C3/C2` is ordered by `(T, C_n)` labels ([F01] eq
  53; [CFP00] eq 46); `T ⊃ D2` uses `(x, y, z)`.
- **Dihedral chain basis = canonical basis up to the swap.** Our canonical
  axial matrices are `D(a) = diag(w^k, w^{−k})`, `D(b) = [[0,1],[1,0]]`
  (ordinary) and `D(s(a)) = diag(u^k, u^{−k})`, `D(s(b)) = [[0,i],[i,0]]`
  (projective) with `a = C_n`, `b = C2'(1)`. Chen's `D_n ⊃ C_n` matrices
  (eq 2.12a) are `diag(e^{−imφ_1}, e^{imφ_1}) = diag(w^{−m}, w^{m})` and
  `(−1)^ξ [[0,1],[1,0]]` for `C2x`: the **same matrices with the two
  components interchanged** (the swap leaves `D(b)` invariant). No phase
  differs. For `C_nv`, `D_nd`, `D_nh` the isomorphisms of [CFMB00] eqs 3.2,
  4.7, 5.4 realise the same matrices on the image elements, which is our
  "type III groups are `R` with parities attached" rule, including `σ(k) ↔
  C2'(k)` for `C_nv` with the plane `σ_x` (normal along `x`) as `C2x`; the
  book's `D(C2x)` "depends on both `m` and `n`" ([CFMB00] eq 7.1 and §VII.A)
  — the same observation that led to the spin-harmonic rule of version 1.
- **Polyhedral chain bases are new bases.** Our canonical polyhedral bases
  diagonalise the polyhedral principal element (`C3^1(1)` for `T`, `C4^1(1)`
  for `O`, `C5^1(1)` for `I`) and freeze the remaining phases from a product
  construction (`E1/2 ⊗ T`, `T1 ⊗ T1`, …). They coincide with the chains
  `T ⊃ C3`, `O ⊃ C4`, `I ⊃ C5` as far as block structure goes (both
  diagonalise the same element), but the phases differ, and for `I5/2` the
  canonical basis (product basis `E1/2 ⊗ T2` ordered by weight `5, 3, 1, −1,
  −3, −5`) is **not** the `I ⊃ D5 ⊃ C5` basis, which separates the two `μ =
  5/2` vectors by the `D5` irreps `[5/2±]` ([FCD99b] eq 37, 44). A `--chain
  I ⊃ D5 ⊃ C5` request must therefore produce a basis different from
  `canonical` for `I5/2` alone; for all other irreps of `I` the `I ⊃ C5` and
  `SO3 ⊃ I` bases coincide ([FCD99b] p. 1057). Likewise `T ⊃ D2`, `T ⊃ C2`,
  `O ⊃ T ⊃ …`, `O ⊃ D4 ⊃ …`, `O ⊃ D3 ⊃ …`, `D_n ⊃ C2x` are bases in which the
  principal element is *not* diagonal.
- **Fields.** The chain bases need larger conductors than the canonical ones
  (`12` for `2T`, `8` for `2O`, `20` for `2I`): `T ⊃ C3` has `√(1/3),
  √(2/3)` (`Q(ζ_24)`); `T ⊃ C2` has `√i/√2` and `ε` (`Q(ζ_24)`); `O ⊃ C4`
  has `√2, √3, i` (`Q(ζ_24)`); `I ⊃ C5` has `c_k, √5, √2, √6` (`Q(ζ_120)`);
  `I ⊃ D5 ⊃ C5` has `s_k, t_k, √8` (`Q(ζ_40)`, `Q(ζ_120)` with `H`). These are
  the conductors of the Altmann & Herzig polyhedral data (24, 24, 120), so
  one cyclotomic field per document remains possible. `D_n ⊃ C_n` needs
  `lcm(2n, 4)` as the canonical basis does; `D_n ⊃ C2x` the same (cosines
  and sines of `2kmπ/n`). The program output form `sqrt(p/q) exp(i cos⁻¹
  (m/n))` ([CP99] abstract) is a printing form for unit-modulus entries such
  as `(1 + i√2)/√3`, which still lie in a cyclotomic field; no entry outside
  `Q(ζ_120)` was seen.
- **CG coefficient fields.** `D_n` CG coefficients add `1/√2` to the matrix
  field (`Q(ζ_{lcm(2n, 8)})`, e.g. `Q(ζ_24)` for `D3` whose matrices live in
  `Q(ζ_12)`); `T ⊃ C3` CG coefficients lie in `Q(ζ_24)` (`1/√3, √(2/3),
  1/√2, 1/√6`); `O ⊃ C4` multiplicity-free coefficients lie in `Q(ζ_24)`, but
  the orthogonalised second copies of `F0 × G3/2 → G3/2` introduce `√5, √10,
  √15` (`Q(ζ_120)`), a consequence of the Gram–Schmidt step, not of the
  group.

### 5.3 Irrep symbols and element labels

- **Modified Mulliken notation** ([CF98c] §III; [CFMB00] §II.A): letter by
  dimension `A, E, F, G, H, I` (`F` for three-dimensional irreps, where we
  use `T`), subscript `m_m = m_max`, the largest weight of the cyclic subgroup
  of the natural chain, except `E−3/2` of `T` (`m_min`). Thus `T`: `A0, A1,
  A−1, F, E1/2, E3/2, E−3/2` (ours `A, 2E, 1E, T, E1/2, 2E3/2, 1E3/2`);
  `O` in [FC99]/[F03]: `A0, A2, E, F0, F2, E1/2, E3/2, G3/2` (ours `A1, A2,
  E, T1, T2, E1/2, E5/2, G3/2`) and in [CFP00]/[F01]: `A1, A2` or `A_σ`,
  `F_σ`, `σ = 1, 2`; `I`: `A, F1, F2, G, H` and d.v. `E1/2, E3/2, G3/2,
  I5/2` where **`E3/2` is our `E7/2`** (A&H `E7/2`, Butler `1/2~`)
  ([FCD99b] Table 1 footnote); `D_n`: `A_{ρν}` / `A_{m0 ν}` (`A_{11} = A1`,
  `A_{12} = A2`, `A_{21} = B1`, `A_{22} = B2` for even `n`; `A_{n/2, ∓i}` the
  one-dimensional d.v. pair for odd `n`) and `E_m` with `m = 1/2, 1, 3/2, …`
  (ours `E_k`, `E_{k/2}`; identical subscripts). [CFMB00] Table I maps these
  to Koster for `D2–D6` and shows Koster's `Ē2, Ē3` of `D6` are `E5/2, E3/2`
  in that order. The programs print Koster `Γ_i` for the irreps and
  `[Γ_j γ_k, …]` for the components ([CP99] §8; [PC99] p. 77).
- **Element labels.** Bradley & Cracknell throughout the cubic papers
  (`C31⁺` … `C34⁺` toward vertices 1–4 of [CF98a] Fig 1, `C2a … C2f`,
  `σ_da … σ_df`, `S6j^±`, `S4x^±`), with the caution that B–C's Table 6.2 has
  a misprint (columns 6 and 8 of row 19 interchanged, [CP99] §8; [PC99]
  p. 68) and that B–C's `D^E(S4^±) = D^E(C4^±)` contradicts their own
  characters `χ(S4^±) = χ(C4^∓)` ([CFMB00] §VII.B). Our canonical labels
  (`C3^1(k)`, `C2'(k)`, `sigma_d(k)`, …) need a map to B–C through the
  element tables ([CP99] Table 3) — a transcription task for #48, akin to the
  A&H label convention.
- **Coordinates.** The papers' frames (`C2` axes of `T`, `O` along `x, y, z`;
  `C5,6` of `I` along `z`) only matter for SAFs; the chain bases themselves
  are coordinate-free, as our abstract-group model requires.

### 5.4 What the papers verify against

- `T ⊃ C3`: differs from Damhus, Harnung & Schäffer (1984) "in phases"
  ([CF98a] p. 488). `T ⊃ D2` SAFs identical to Bradley & Cracknell pp. 66–70
  ([CF98b] p. 5509). `T ⊃ C2` SAFs agree with Damhus for `j = 1/2, 3/2`
  ([CF98b] p. 5515).
- `O ⊃ C4`: identical to Butler pp. 527–531 for `j ≤ 8` in multiplicity-free
  cases, differ with multiplicity ([FC99] p. 20). `O ⊃ D4 ⊃ D2`: identical to
  B–C up to overall phase, multiplicity cases differ except `F1, F2`; irrep
  matrices identical to B–C ([CFP00] p. 597). `O ⊃ T ⊃ D2` CG and irreps of
  the program identical to Altmann & Herzig ([CP99] §9).
- `I ⊃ C5`: identical to Damhus p. 439 for `l = 0–5`, differ from Butler by
  phases (coordinate axes), real unlike McLellan ([FCD99a] §9); `I ⊃ D5 ⊃
  C5`: identical to Damhus for `j = 1/2–9/2` multiplicity-free, identical to
  Butler up to phases, all real except `I5/2` at `±5/2` ([FCD99b] §8).
- Dihedral: identical to Altmann & Herzig for one-dimensional irreps, phase
  differences for two-dimensional ones; A&H "have not chosen the irreducible
  matrices and CG coefficients of isomorphic groups identical" ([CFMB00]
  §VII.A); the factor system of `D6†` built from the `E1/2` matrices equals
  Altmann's p. 253 ([CFMB00] §VII.B).

---

## Summary for the map

The eigenfunction method is the simultaneous diagonalisation of a class
operator combination `C` (CSCO-I), the chain operator `C(s) = aC_H + bC_K`
(CSCO-II) and its intrinsic copy `C̄(s)` (CSCO-III) in the group algebra,
whose eigenvectors are the chain-adapted irrep matrix elements; CG
coefficients are the eigenvectors of the CSCO-II in the Kronecker product.
The double-induced technique reduces `G ⊃ C_n` to a `q × q` problem over the
double-coset generators (`q = 2, 2, 3, 4` for `D_n, T, O, I`), the
factorization lemma solves `T ⊃ D2` and all dihedral groups without a group
table, and the right-induced technique builds `O ⊃ T ⊃ …` from `T`. Printed
exact matrices exist for `T ⊃ C3`, `T ⊃ D2`, `T ⊃ C2`, `O ⊃ C4`, `O ⊃ T ⊃
C3`, `O ⊃ T ⊃ D2 ⊃ C2`, `I ⊃ C5`, `I ⊃ D5 ⊃ C5` and (algebraically, any `n`)
`D_n ⊃ C_n`, `D_n ⊃ C2x` with `C_nv`, `D_nd`, `D_nh` by isomorphism; printed
CG coefficients for `T ⊃ C3`, `T ⊃ D2 ⊃ C2`, `T ⊃ C2`, `O ⊃ C4`, `O ⊃ T ⊃ C3`
(d.v. × d.v. only) and `D_n ⊃ C_n` (algebraic). The eigenproblem leaves one
phase per basis vector and the multiplicity basis open; the papers close them
by external conventions (agreement with Butler, real SAFs, hermiticity,
associated-irrep relations), while the CG phase is fixed by a positive
coefficient of a chosen trial product state and the multiplicity by the
trial state's intrinsic quantum numbers (copies not orthogonal). Against
version 1: `m = −e` (components reversed, Chen's `E3/2` of `T`, `O`, `I` are
our `2E3/2`, `E5/2`, `E7/2`), the dihedral untilded set is the `[0, 2π)`
section, the `O†` table is Herzig's with `C2b, C2c, C2f` reversed, `F` names
three-dimensional irreps, labels are Bradley–Cracknell or Koster, and the
chain bases need conductors 24 (`T`, `O`) and 120 (`I`), with `D_n` CG
coefficients needing `lcm(2n, 8)` and orthogonalised multiplicities in `O`
needing `√5`.
