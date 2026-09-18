# Irrep labels

Decided in wayfinder ticket
[Irrep labelling rules: Mulliken symbols for arbitrary n and Bethe aliases](https://github.com/hydeik/pointgroup-tables/issues/7).
Builds on [Abstract groups](./abstract-groups.md) (every element is a pair
`(r, parity)` with `r` in a rotation group `R`),
[Double groups](./double-groups.md) (the canonical section `s` and the
projective irreps) and [Label conventions](./label-conventions.md) (the
two-layer labelling model, reused here for irreps).

## Two layers

Every irrep carries a **canonical Mulliken symbol**, derived by the fixed
rules below from its dimension and from its characters on a few named
elements of the abstract group. It is convention-free, defined for every
axis order `n`, ordinary and projective, and always present.

A **reference irrep label** is the spelling a reference work (Altmann &
Herzig, ...) uses for the same irrep. It is produced by the selected
**label convention**, from rewrite rules on the canonical symbol plus
per-group overrides, exactly as for element labels. Reference irrep labels
appear only when a convention is selected.

A **Bethe symbol** `Γ_i` is a third name, an alias present only for the 32
crystallographic point groups; see [Bethe aliases](#bethe-aliases).

Irrep names are **output only**. Tables reference irreps by their index in
the canonical irrep order, never by name.

## Canonical symbol grammar

ASCII only, no whitespace:

```
symbol    := component? letter subscript? suffix?
component := "1" | "2"                      # member of a complex-conjugate pair
letter    := "A" | "B" | "E" | "T" | "G" | "H" | "I"
subscript := digit+                         # secondary subscript (A, B, T) or k (E)
           | odd "/2"                       # projective subscript j = k/2
suffix    := "g" | "u" | "'" | "''"
```

Regular expression accepting exactly these symbols:

```
^[12]?[ABETGHI](\d+|\d*[13579]/2)?(g|u|'|'')?$
```

Examples: `A`, `A1g`, `B3u`, `E`, `E2`, `1E1`, `2E`, `E1/2`, `1E3/2u`,
`A1'`, `E''`, `T2g`, `G3/2`, `I5/2g`. The letter, subscript and suffix are
also emitted as separate fields of the [structured irrep record](#structured-irrep-record),
so no consumer needs to parse the symbol.

## Named elements

The rules refer to a handful of elements of the abstract group, all fixed by
the [canonical element order](./abstract-groups.md) and independent of any
coordinate system.

| Name | Element | Meaning by family |
| ---- | ------- | ----------------- |
| **principal element** `c` | `(1, 0)` of `R`, with its parity in `G`; order `N` | `C_n` for `Cn`, `Cnv`, `Dn`, `Cnh`/`Dnh` (even `n`), `S2n`/`Dnd` (odd `n`); `S_2n` for `S2n`/`Dnd` (even `n`); `S_n^-1` (order `2n`) for `Cnh`/`Dnh` (odd `n`) |
| **proper principal element** `c₀` | `c` itself, except for `Cnh`/`Dnh` with odd `n`, where `c₀ = (2, 0) = C_n` (order `n`) | the principal rotation of `G⁺` |
| **secondary element** `s` | `(0, 1)` of `R = D_N`, with its parity | `C2'(1)` in `Dn`, `Dnh`, `Dnd`; `sigma_v(1)` in `Cnv` |
| **inversion** `i` | `(e, −1)`, type II only | |
| **horizontal reflection** `σh` | `(N/2, −1)` of `R = C_2n` or `((N/2, 0), −1)` of `R = D_2n`, `Cnh`/`Dnh` with odd `n` only | `sigma_h` |
| **polyhedral principal element** `c` | `T`, `Th`: `C3^1(1)`, the positive generator of the first `C3` axis; `O`, `Oh`, `Td`: the order-4 element of `R` with the smallest canonical index, with its parity (`C4^1(1)`, or `S4^-1(1)` in `Td`); `I`, `Ih`: `C5^1(1)`, the positive generator of the first `C5` axis. In every case the element of the principal order (`3`, `4`, `5`) with the smallest canonical index | |

`N` is the order of `c`, `ζ_N = exp(2πi/N)`, and `d` is the dimension of the
irrep. For a projective irrep the characters are those of the lifted
elements `s(c)`, `s(s)`, ... in the canonical section.

## Rules: axial families

Let `Cn`, `Cnv`, `Cnh`, `S2n`, `Dn`, `Dnd`, `Dnh` be the axial families with
their rotation group `R` from the [family table](./abstract-groups.md).

### Letter and `k` (ordinary irreps)

1. **One-dimensional, real on `c₀`.** `χ(c₀) = +1` gives `A`, `χ(c₀) = −1`
   gives `B`. `B` occurs only when `N` is even.
2. **One-dimensional, complex on `c₀`.** `χ(c₀) = ζ_N^{±k}` with
   `1 ≤ k < N/2`. The irrep and its complex conjugate form a pair; both are
   `E_k`, distinguished by the component prefix (rule 8).
3. **Two-dimensional.** The eigenvalues of `c₀` are `ζ_N^{±k}` with
   `1 ≤ k < N/2`; the irrep is `E_k`.
4. **Subscript elision.** When only `k = 1` is possible (`N ≤ 4`) the
   subscript is omitted: `E`, `1E`, `2E`. For `N ≥ 5` every `E` carries its
   `k` (`E1`, `E2`, ...).

For `Cnh`/`Dnh` with odd `n` the rules use `c₀ = C_n` (order `n`), so their
ordinary irreps are named as in `Cn`/`Dn` and then primed by rule 7.

### Secondary subscript

5. **`1` / `2`.** In the families with a secondary element (`Cnv`, `Dn`,
   `Dnd`, `Dnh`), an `A` or `B` irrep gets subscript `1` when `χ(s) = +1`
   and `2` when `χ(s) = −1`. `E` irreps carry `k`, never a secondary
   subscript. Consequences: in `D4`, `B1` is symmetric under `C2'` and `B2`
   under `C2''`; in `C4v`, `B1` is symmetric under `sigma_v` and `B2` under
   `sigma_d`; in `D2d`, `B1` is symmetric under `C2'` and `B2` under
   `sigma_d`.
6. **`D2` and `D2h`.** `D2` has no unique principal axis and is the one
   exception: `B1` is the irrep with `χ(C2) = +1` (the principal element
   `(1, 0)`), `B3` the one with `χ(C2'(1)) = +1`, `B2` the one with
   `χ(C2''(1)) = +1`. This is the standard assignment `B1 ~ z`, `B2 ~ y`,
   `B3 ~ x` once `C2'(1)` is taken along `x`, as in the `D4`, `D6`
   conventions of Altmann & Herzig. `D2h` adds `g`/`u`.

### Suffix

7. **`g` / `u`** (type II groups): `χ(i) = +d` gives `g`, `−d` gives `u`.
   **`'` / `''`** (`Cnh`, `Dnh` with odd `n`): `χ(σh) = +d` gives `'`,
   `−d` gives `''`. Projective irreps of `Cnh`/`Dnh` with odd `n` carry no
   prime (their `j` already encodes the behaviour under `σh`); projective
   irreps of type II groups do carry `g`/`u`.

### Complex-conjugate pairs

8. **Component prefix.** The two members of a complex-conjugate pair of
   one-dimensional irreps (or, in `T`, `Th` and their double groups, of
   two-dimensional irreps) are separate irreps with the same letter,
   subscript and suffix, prefixed `1` and `2`. Member `1` is the one whose
   character at the first element, in canonical order, where the two
   characters differ has positive imaginary part. For the axial families
   this is the member with the **positive exponent**: `χ(c₀) = ζ_N^{+k}`
   (ordinary), `χ(s(c)) = ζ_2N^{+k}` (projective, `k < N`), or
   `χ(s(s)) = +i` (the projective pair with `k = N`, odd `N`, in `Dn`,
   `Cnv`, `Dnd`, `Dnh`). Each member records the other as its `partner`.
   No merged "real E" row is ever emitted.

### Projective irreps

9. **Subscript `j = k/2`.** The eigenvalues of `s(c)` (the lift of the
   principal element `c`, of order `N`; note `c`, not `c₀`) are
   `ζ_2N^{±k}` with `k` odd and `1 ≤ k ≤ N`. The irrep is `E_{k/2}`,
   written `Ek/2` (`E1/2`, `E3/2`, `1E5/2`, ...). The subscript is never
   elided.
10. **Dimension and pairing.** `k < N` gives, in `Cn`-based groups
    (`R = C_N`), a complex pair `1Ek/2`, `2Ek/2`, and in `Dn`-based groups
    (`R = D_N`) a two-dimensional `Ek/2`. `k = N` (only when `N` is odd)
    gives one-dimensional irreps with `χ(s(c)) = −1`: in `Dn`-based groups
    the pair `1EN/2`, `2EN/2` (`χ(s(s)) = ±i`); in `Cn`-based groups a
    single real irrep, named by rule 11.
11. **Real one-dimensional projective irreps** (`C1`, `Ci`, and `Cn`, `S2n`
    with `R = C_N`, `N` odd): rule 1 applied to `s(c)`. `C1` and `Ci` have
    `χ(s(E)) = +1` and get `A1/2` (`A1/2g`, `A1/2u`); for odd `N ≥ 3` the
    irrep with `k = N` has `χ(s(c)) = −1` and gets `BN/2`: `B3/2` in `C3`,
    `B3/2g`, `B3/2u` in `S6`, `B5/2` in `C5`. This is the `B` of
    `Cn† ≅ Z_2N` read as `C_2N`. Bradley & Cracknell (and the Bilbao tables,
    which follow them) write these as `Ā` and `Ē`; whatever Altmann & Herzig
    write is recorded by the transcription task and mapped by the label
    convention.

Rule 9 applied with `c = S_n^-1` (order `2n`) in `Dnh`/`Cnh` with odd `n`
gives `D3h†`: `E1/2` (`χ(S3) = √3`), `E3/2` (`0`), `E5/2` (`−√3`); with
`c = S_2n` in `D2d`: `E1/2` (`χ(S4) = √2`), `E3/2` (`−√2`); these agree with
Altmann & Herzig.

## Rules: polyhedral families

Dimension fixes the letter of every irrep, ordinary or projective:
`1 → A`/`B`, `2 → E`, `3 → T`, `4 → G`, `5 → H`, `6 → I`. `F` is never a
canonical letter: Altmann & Herzig and Bradley & Cracknell both write the
four-dimensional projective irreps as `F` (`F3/2`, `F̄`) and the
two-dimensional projective pair of `T` as `1F3/2`, `2F3/2`; the canonical
symbols keep the dimension rule and the label convention restores `F`.

12. **`T`, `Th`.** `A`, the complex pair `1E`, `2E` (`χ(C3^1(1)) = ζ_3^{±1}`,
    member `1` positive exponent), `T`; no subscripts. `Th` adds `g`/`u`.
13. **`O`, `Oh`, `Td`.** `A1`, `A2`, `E`, `T1`, `T2`; the subscript of `A`
    and `T` is `1` when `χ(c) = +1` on the polyhedral principal element
    (`C4^1(1)` in `O`/`Oh`, `S4^-1(1)` in `Td`) and `2` when `χ(c) = −1`.
    `Oh` adds `g`/`u`.
14. **`I`, `Ih`.** `A`, `T1`, `T2`, `G`, `H`; `T1` is the three-dimensional
    irrep with `χ(C5^1(1)) = φ = (1 + √5)/2` and `T2` the one with
    `χ(C5^1(1)) = 1 − φ` (exact comparison in `Q(√5)`). `Ih` adds `g`/`u`.
    Because the outer automorphism of `I` swaps the two classes of
    five-fold rotations, this rule and the canonical label `C5^1(1)` fix
    each other; a reference whose `C5` is our `C5^2` swaps `T1`/`T2` through
    its label convention.
15. **Projective irreps.** The letter follows the dimension; the subscript
    `j` is the smallest half-integer such that the irrep occurs in the
    restriction of the SU(2) irrep `D^j` to `G†`, which the characters below
    fix exactly (all signs are those of the canonical section, under which
    every lifted rotation has `w > 0`):

    | Group | Irrep | Fixing character |
    | ----- | ----- | ---------------- |
    | `T†`, `Th†` | `E1/2` | two-dimensional, self-conjugate, `χ(s(C3^1(1))) = 1` |
    | | `1E3/2`, `2E3/2` | two-dimensional pair, `χ(s(C3^1(1))) = ζ_3^{±1}`, member `1` positive exponent |
    | `O†`, `Oh†`, `Td†` | `E1/2` | two-dimensional, `χ(s(c)) = +√2` |
    | | `E5/2` | two-dimensional, `χ(s(c)) = −√2` |
    | | `G3/2` | four-dimensional |
    | `I†`, `Ih†` | `E1/2` | two-dimensional, `χ(s(C5^1(1))) = φ` |
    | | `E7/2` | two-dimensional, `χ(s(C5^1(1))) = 1 − φ` |
    | | `G3/2` | four-dimensional |
    | | `I5/2` | six-dimensional |

    Type II polyhedral double groups add `g`/`u` by rule 7. The outer
    automorphism of `I` swaps `E1/2` with `E7/2` exactly as it swaps `T1`
    with `T2`, so rule 14 and this rule are tied to the same class
    `C5^1(1)`.

## Canonical irrep order

Every table lists irreps in this order, which depends on nothing but the
canonical symbols:

1. ordinary irreps first, then projective irreps;
2. within each block, by suffix: none, or `g` before `u`, or `'` before
   `''`;
3. within a suffix block, by letter in the order `A`, `B`, `E`, `T`, `G`,
   `H`, `I`, then by subscript ascending (`1`, `2`, `3`; `k`; `k/2`), then
   component `1` before `2`.

`D4h`: `A1g A2g B1g B2g Eg A1u A2u B1u B2u Eu`, then `E1/2g E3/2g E1/2u
E3/2u`. `C6`: `A B 1E1 2E1 1E2 2E2`, then `1E1/2 2E1/2 1E3/2 2E3/2 1E5/2
2E5/2`. `D3h`: `A1' A2' E' A1'' A2'' E''`, then `E1/2 E3/2 E5/2`. `O`:
`A1 A2 E T1 T2`, then `E1/2 E5/2 G3/2` (the Altmann & Herzig order). `I`:
`A T1 T2 G H`, then `E1/2 E7/2 G3/2 I5/2`. The Bethe number is not an
ordering key: Koster's `Γ_i` do not follow this order (`D3h`: `E'` is `Γ6`,
`E''` is `Γ5`).

## Examples

| Group | Ordinary irreps | Projective irreps |
| ----- | --------------- | ----------------- |
| `C1` | `A` | `A1/2` |
| `Cs` | `A' A''` | `1E1/2 2E1/2` |
| `C2` | `A B` | `1E1/2 2E1/2` |
| `C3` | `A 1E 2E` | `B3/2 1E1/2 2E1/2` |
| `C4` | `A B 1E 2E` | `1E1/2 2E1/2 1E3/2 2E3/2` |
| `S4` | `A B 1E 2E` | `1E1/2 2E1/2 1E3/2 2E3/2` |
| `S6` | `Ag 1Eg 2Eg Au 1Eu 2Eu` | `B3/2g 1E1/2g 2E1/2g B3/2u 1E1/2u 2E1/2u` |
| `C3h` | `A' 1E' 2E' A'' 1E'' 2E''` | `1E1/2 2E1/2 1E3/2 2E3/2 1E5/2 2E5/2` |
| `C5` | `A 1E1 2E1 1E2 2E2` | `B5/2 1E1/2 2E1/2 1E3/2 2E3/2` |
| `D2` | `A B1 B2 B3` | `E1/2` |
| `C3v`, `D3` | `A1 A2 E` | `E1/2 1E3/2 2E3/2` |
| `C4v`, `D4`, `D2d` | `A1 A2 B1 B2 E` | `E1/2 E3/2` |
| `D3h` | `A1' A2' E' A1'' A2'' E''` | `E1/2 E3/2 E5/2` |
| `D5` | `A1 A2 E1 E2` | `E1/2 E3/2 1E5/2 2E5/2` |
| `D6h` | `A1g A2g B1g B2g E1g E2g A1u ... E2u` | `E1/2g E3/2g E5/2g E1/2u E3/2u E5/2u` |
| `D4d` | `A1 A2 B1 B2 E1 E2 E3` | `E1/2 E3/2 E5/2 E7/2` |
| `T` | `A 1E 2E T` | `E1/2 1E3/2 2E3/2` |
| `O`, `Td` | `A1 A2 E T1 T2` | `E1/2 E5/2 G3/2` |
| `Oh` | `A1g A2g Eg T1g T2g A1u A2u Eu T1u T2u` | `E1/2g E5/2g G3/2g E1/2u E5/2u G3/2u` |
| `I` | `A T1 T2 G H` | `E1/2 E7/2 G3/2 I5/2` |

## Structured irrep record

Every irrep is emitted with its symbol and a decomposition:

| Field | Values |
| ----- | ------ |
| `index` | position in the canonical irrep order |
| `symbol` | canonical Mulliken symbol |
| `reference_symbol` | present only when a label convention is selected |
| `bethe` | the Bethe symbol as the ASCII string `Gamma<i>` or `Gamma<i>+` / `Gamma<i>-` (`"Gamma4"`, `"Gamma6-"`); present only for the 32 crystallographic groups |
| `dimension` | `d` |
| `projective` | `true` or `false` |
| `letter` | `A`, `B`, `E`, `T`, `G`, `H`, `I` |
| `component` | `1` or `2`; absent unless the irrep is a member of a complex-conjugate pair |
| `partner` | index of the other member of the pair; absent otherwise |
| `subscript` | the secondary subscript or `k` as an integer, or `k/2` as the string `"3/2"`; absent when the symbol has none |
| `suffix` | `g`, `u`, `'`, `''`; absent when the symbol has none |
| `frobenius_schur` | `1`, `0`, `-1` |

## Bethe aliases

The Bethe symbol `Γ_i` of an irrep follows Koster, Dimmock, Wheeler and
Statz, *Properties of the Thirty-Two Point Groups* (1963), which numbers the
ordinary irreps and then the projective irreps of each crystallographic
point group, with a `±` superscript for the parity under inversion in the
centrosymmetric groups. No rule generates these numbers for arbitrary `n`,
so they are a transcribed table, keyed by group symbol and canonical
Mulliken symbol, shipped with the package
(`src/pointgroup_tables/data/conventions/bethe/koster.toml`, generated from
the transcription by `tools/koster_bethe.py generate`; each entry carries the
irrep's characters on the group's identifying elements, by which the tool
attaches it to the computed irrep) and emitted as the `bethe` field
(spelled `Gamma4`, `Gamma6-`) for the 32 crystallographic groups and their
double groups whatever label convention is selected. Groups outside the 32
(any other `n`, `I`, `Ih`) have no `bethe` field. Bradley & Cracknell
(Table 2.2) and the Bilbao point-group tables (which print the Koster
numbers) serve as second sources when the transcription is checked; the
Bilbao single-group and double-group pages disagree with each other for
`C6`, so Koster is the authority. The transcription is the task
[Transcribe Koster Γ_i numbering](https://github.com/hydeik/pointgroup-tables/issues/15).

## Reference irrep labels in a label convention

The label convention file gains an `[irreps]` table:

```toml
[irreps]
component_style = "prefix"          # 1E stays 1E ; "superscript" gives ^1E
fraction_style = "plain"            # E1/2 stays E1/2 ; "subscript" gives E_{1/2}
projective_parity_separator = ","   # E1/2g -> E1/2,g ; "" keeps E1/2g (ordinary A1g is never separated)
letters = { G = "F" }               # ordinary irreps: canonical -> reference letter
projective_letters = { G = "F", B = "A" }   # projective irreps only: G3/2 -> F3/2, B3/2 -> A3/2

[irreps.overrides.T]                # per-group explicit map, canonical -> reference
"1E3/2" = "1F3/2"
"2E3/2" = "2F3/2"
```

Rules rewrite the spelling; `letters` and `projective_letters` rename a
letter for every group; overrides handle a reference that swaps the members
of a pair, keeps a parent letter for a split pair, or names an irrep
irregularly. The Altmann & Herzig table above is the one the transcription
of the irrep symbol rows fixed
([Altmann & Herzig label-convention rules and overrides](https://github.com/hydeik/pointgroup-tables/issues/14)):
the book writes every four-dimensional irrep as `F` (`F` of `I`, `F3/2`,
`F3/2,g`), the real one-dimensional projective irreps of `C3` and `S6` as
`A3/2` (`A3/2,g`, `A3/2,u`), the projective `g`/`u` after a comma
(`E1/2,g`, `1E1/2,u`), and the pair of `T` / `Th` as `1F3/2`, `2F3/2`
(`1F3/2,g`, ... in `Th`); `A1/2`, `A'`, `1E1/2`, `E1/2`, `I5/2` and every
ordinary symbol are spelled as the canonical rules do. Whether the book
assigns `1E` / `2E`, `B1` / `B2` and `T1` / `T2` as the canonical rules do
is a matter of characters, not spelling: it is checked by `REF-5` against
the transcribed character tables
([Transcribe Altmann & Herzig character tables and sample axial irrep matrices](https://github.com/hydeik/pointgroup-tables/issues/19),
extended to every axial group with `n ≤ 10` by
[Reference irrep symbols of the `altmann-herzig` convention beyond the transcribed character tables](https://github.com/hydeik/pointgroup-tables/issues/36)),
and any disagreement becomes an `[irreps.overrides]` entry. `C2v` and `C6v`
were settled ahead of that by the element isomorphism of the label
convention ([Label conventions](./label-conventions.md#the-altmann--herzig-convention)).

## Verification

An implementation must check, for every group it constructs:

- every irrep receives exactly one canonical symbol matching the grammar,
  and the symbols of a group are pairwise distinct;
- the symbol is a function of the character only (recomputing it from the
  character table reproduces it);
- `dimension` equals the character of the identity, `frobenius_schur` is
  `(1/|G|) Σ_g ω(g, g) χ(g²)` (with `ω ≡ 1` for ordinary irreps), and a
  `partner` is present exactly when `frobenius_schur = 0` and the
  conjugate irrep is a different irrep;
- against Altmann & Herzig, for the groups transcribed by the transcription
  task: the canonical symbol, after the Altmann & Herzig convention, equals
  the book's symbol row by row.

Comparison of the characters themselves belongs to the verification
chapter.
