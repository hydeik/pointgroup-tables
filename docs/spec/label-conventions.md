# Label conventions

Decided in wayfinder ticket
[Label conventions: mapping abstract elements to geometric labels](https://github.com/hydeik/pointgroup-tables/issues/3).
Builds on [Abstract groups](./abstract-groups.md): every element is a pair
`(r, parity)` with `r` in a rotation group `R`. Prototype of the canonical
rules: `prototypes/abstract_group_prototype.py labels <symbols>` on branch
`prototype/abstract-group`.

## Two layers

Every element carries a **canonical label**, derived by a fixed formula from
`(r, parity)`, the family, the group type and `H`. It is convention-free,
defined for every axis order `n`, and always present.

A **reference label** is the name a reference (Altmann & Herzig, ...) gives
the same element. It is produced by a named **label convention**: a data file
of rewrite rules applied to the canonical label plus, where the reference is
irregular, a per-group override table. Reference labels appear only when a
convention is selected; at most one convention is selected per run.

Labels are **output only**. Tables reference elements by their index in the
canonical element order, never by label.

## Canonical label grammar

ASCII only. `m ≥ 3`, `j` a signed integer with `1 ≤ |j| ≤ m/2` and
`gcd(j, m) = 1`, `k` a 1-based index. A label never contains whitespace.

```
label     := "~"? base
base      := "E" | "i" | rotation | secondary | improper | plane
rotation  := "C2" | "C" m "^" j | "C" m "^" j "(" k ")" | "C2(" k ")"
secondary := "C2'(" k ")" | "C2''(" k ")"
improper  := "S" m "^" j | "S" m "^" j "(" k ")"
plane     := "sigma_h" | "sigma_h(" k ")" | "sigma_v(" k ")"
           | "sigma_d(" k ")" | "sigma(" k ")"
```

Regular expression accepting exactly these labels:

```
^~?(E|i|C2|C2\(\d+\)|C2''?\(\d+\)|C\d+\^-?\d+(\(\d+\))?|S\d+\^-?\d+(\(\d+\))?|sigma_h(\(\d+\))?|sigma_[vd]\(\d+\)|sigma\(\d+\))$
```

The optional `~` prefix marks a barred element of the double group (`~E` is
the full turn Ē); see [Double groups](./double-groups.md). It appears only in
double-group tables.

## Canonical rules: axial families

Let `N` be the order parameter of `R` (`C_N` or `D_N`), `n′` the axis order of
the proper subgroup `G⁺` (`n′ = n`, except `S2n` where `n′ = n`; for type III
groups with `R = C_2n`/`D_2n` the proper elements are the even `k`, and
`k′ = k / 2` indexes them within `G⁺`). Turns are fractions of a full
rotation.

1. **Rotation** `(k, 0, +1)`: `k = 0` gives `E`. Otherwise reduce `k/N` to
   `j/m` in lowest terms, replace `j` by `j − m` when `2j > m` so that
   `|j| ≤ m/2`; label `C{m}^{j}`, or `C2` when `m = 2`.
2. **Secondary axis** `(k, 1, +1)` with `k′` as above: `n′` odd gives
   `C2'(k′+1)`; `n′` even gives `C2'(k′/2 + 1)` for even `k′` and
   `C2''((k′−1)/2 + 1)` for odd `k′`. Even `k′` is *defined* to be the primed
   class.
3. **Improper rotation** `(k, 0, −1)`: this is `i · C(k/N)`, i.e. a
   rotoreflection through the turn `g = k/N + 1/2 (mod 1)`. `g = 0` gives
   `sigma_h`, `g = 1/2` gives `i`, otherwise reduce `g` as in rule 1 and label
   `S{m}^{j}`. Hence `i·C3 = S6^-1`, never `S6^5`.
4. **Vertical plane** `(k, 1, −1)`, the plane perpendicular to secondary
   direction `k`:
   - `Cnv`: `n` odd gives `sigma_v(k+1)`; `n` even gives `sigma_v(k/2 + 1)`
     for even `k` and `sigma_d((k−1)/2 + 1)` for odd `k`;
   - `Dnh`: the plane contains secondary direction `d = k + N/2 (mod N)`;
     with `d′ = d / step` (the index within `G⁺`) the plane is `sigma_v(x)`
     if direction `d′` is a `C2'(x)` axis and `sigma_d(x)` if it is a
     `C2''(x)` axis (rule 2). In `D6h` the planes `i·C2'` are therefore
     `sigma_d`, and `sigma_v(x)` always contains `C2'(x)`;
   - `Dnd`: every vertical plane is `sigma_d(x)`, indexed by `k + 1` (`n`
     odd, `R = D_n`) or `(k−1)/2 + 1` (`n` even, `R = D_2n`, odd `k`).

Because planes are indexed by the axis they contain, plane labels inside the
improper block do not follow the index order of `R`; the element order itself
is unchanged.

## Canonical rules: polyhedral families

Let `R` be `T`, `O` or `I` with its canonical element order.

5. **Axes.** The axis of `r ≠ E` is the maximal cyclic subgroup containing
   `r`. Axes are numbered per maximal order, in order of first appearance of
   any of their elements in the canonical order. The **positive generator** of
   an axis is its maximal-order element with the smallest canonical index;
   every element of the axis is `gen^j`, giving the turn `j/o`.
6. **Rotations** `(r, +1)`: reduce the turn as in rule 1 and label
   `C{m}^{j}(k)` with `k` the axis number; order-2 elements are `C2(k)`, or
   `C2'(k)` when the element lies on an order-2 axis *and* the group also has
   order-2 elements on higher-order axes (`O`: `3 C2` on the `C4` axes and
   `6 C2'`). `C4^2` on axis `k` is `C2(k)`.
7. **Improper elements** `(r, −1)` (type II: `i·r`; `Td`: `r ∉ T`): `r = E`
   gives `i`; `r` of order 2 gives a reflection named by family: `Th` and
   the `C4`-axis planes of `Oh` are `sigma_h(k)`, the `C2'`-axis planes of `Oh`
   and all planes of `Td` are `sigma_d(k)`, planes of `Ih` are `sigma(k)`;
   otherwise apply rule 3 to the turn and append `(k)`: `i·C3 = S6^-1(k)`,
   `i·C4 = S4^-1(k)`, `i·C5 = S10^-3(k)`, `i·C5^2 = S10^-1(k)`.

Resulting class structure (checked in the prototype): `Td`: `E, 3C2, 8C3,
6sigma_d, 6S4`; `Oh`: `E, 3C2, 6C2', 8C3, 6C4, i, 3sigma_h, 6sigma_d, 8S6,
6S4`; `Ih`: `E, 15C2, 20C3, 12C5, 12C5^2, i, 15sigma, 20S6, 12S10^3, 12S10`.

## Structured element record

Every element is emitted with its label and a decomposition so consumers need
not parse the label:

| Field | Values |
| ----- | ------ |
| `index` | position in the canonical order |
| `label` | canonical label |
| `reference_label` | present only when a label convention is selected |
| `parity` | `1` or `-1` |
| `kind` | `identity`, `rotation`, `inversion`, `reflection`, `rotoreflection` |
| `order` | `m` of the rotation or rotoreflection (`2` for reflections, `1` for `E`, `2` for `i`) |
| `power` | `j`, absent for order ≤ 2 |
| `axis_class` | `principal`, `secondary`, `secondary_double_prime`, `polyhedral`, or absent |
| `axis_index` | `k`, absent when the label has none |
| `inverse` | index of the inverse element |
| `class` | index of the conjugacy class containing the element ([Conjugacy classes](./conjugacy-classes.md#class-records)) |
| `barred` | `true` or `false`; double-group tables only |
| `base_index` | index of the point-group element this element covers; double-group tables only |

Key order and presence rules are fixed in
[Output schema](./output-schema.md#the-elements-table).

## Label conventions (reference labels)

One TOML file per convention, shipped with the package as
`src/pointgroup_tables/data/conventions/labels/<name>.toml` (the namespace
of the [basis conventions](./basis-conventions.md) is `.../conventions/basis/`).
The default is the canonical labelling itself (no file). A convention is
applied to the canonical label of every element by the rules below; a
per-group `overrides` entry wins over every rule.

```toml
name = "altmann-herzig"
source = "S. L. Altmann and P. Herzig, Point-Group Theory Tables, 2nd ed."
exponent_style = "plusminus"      # C3^1 -> C3+, C3^-1 -> C3-, C8^3 -> C8^3+, S10^-3 -> S10^3-
barred_format = "~{label}"        # how the reference spells a barred element; "{label}~" for a suffix

[coverage]                        # the reference's range: verified here, extrapolated beyond
axial_max_n = 10                  # every axial family up to this n
polyhedral = true                 # all seven polyhedral groups

[secondary]                       # the reference's numbering of the secondary directions of D_Nh
"6" = { primed = [0, 4, 8], double_primed = [3, 7, 11] }

[overrides.D2h]                   # per-group explicit map, canonical -> reference
"C2'(1)" = "C2x"

[barred.O]                        # per-group: elements whose reference lift is our barred element
swapped = ["C2'(1)", "C2'(3)", "C2'(5)"]

[irreps]                          # reference irrep symbols, see Irrep labels
letters = { G = "F" }
```

The rule vocabulary of version 1 is exactly `exponent_style`,
`barred_format`, the `[secondary]` table, `[overrides]`, `[barred]`, the
`[irreps]` table of [Irrep labels](./irrep-labels.md#reference-irrep-labels-in-a-label-convention)
and the `[classes.overrides]` table of
[Conjugacy classes](./conjugacy-classes.md#reference-class-labels). The rules
`swap_primes`, `axis_index_shift` and `sigma_spelling` foreseen when this
framework was designed are needed by no group of Altmann & Herzig and are
not part of the vocabulary. A group without an override table still
receives rule-based reference labels.

### Principal elements

`exponent_style = "plusminus"` rewrites `C{m}^{j}` and `S{m}^{j}` as
`C{m}+` / `C{m}-` when `|j| = 1` and as `C{m}^{|j|}+` / `C{m}^{|j|}-`
otherwise; `E`, `i`, `C2`, `sigma_h` and the `~` prefix are unchanged.
`"canonical"` keeps the canonical spelling.

### Secondary elements: the master group and its direction table

Every axial group with `R = D_N` sits inside its **master group** `D_Mh`,
`M = N` for even `N` and `M = 2N` for odd `N`, whose secondary directions
are `0 … M−1` in units of `π/M`; direction `k` of `G` (units `π/N`) is
direction `k·M/N` of the master. A reference numbers the master's
directions once, in `[secondary."M"]`: `primed` lists the poles of its
`C2'1, C2'2, …` and `double_primed` those of `C2''1, C2''2, …`, as integers
in units of `180/M` degrees (`0 ≤ pole < 2M`; the direction is the pole
modulo `M`). The reference label of a secondary element of `G` in direction
`k` is

- for a proper element `(k, 1, +1)`: the master's label of that direction,
  `C2'{i}` or `C2''{i}` with letter and index from the table entry (so the
  four secondary axes of `D4d`, canonical `C2'(1)`, `C2'(2)`, `C2''(1)`,
  `C2''(2)`, are `C2'1 … C2'4` of `D8h`);
- for a plane `(k, 1, −1)`: the canonical letter (`sigma_v` or `sigma_d`)
  with the index `i` of the master's secondary axis in direction `k`, the
  pole of the plane (Altmann & Herzig write every improper operation as
  `i` times the proper operation with the same parameters).

**The twist.** In `Cnv` with `n ≡ 2 (mod 4)` canonical direction `0`
corresponds to the master's direction `M/2` (the direction of `C2''1`),
i.e. direction `k` of `G` is looked up as `k + N/2`. The reference's
`sigma_v` planes of these groups contain the master's `C2'` axes, whereas
the canonical `sigma_v(1)` is perpendicular to direction `0`; the twist
keeps the letters of the two layers equal and, with it, rule 5 of
[Irrep labels](./irrep-labels.md) (`B1` symmetric under `sigma_v`) agrees
with the reference without an irrep override
([ADR-0007](../adr/0007-abstract-element-identification.md)).

**Coverage.** A group whose master `M` has no `[secondary]` entry receives
the spelling rules only and keeps the canonical indices. The convention
file declares its range in `[coverage]` (`axial_max_n`, `polyhedral`), and
the document metadata reports `coverage = "reference"` for a group inside
it and `"extrapolated"` outside ([Output schema](./output-schema.md#metadata)).
The `elements` table is unchanged either way.

### Reference section

For the axial groups the barred assignment of the reference follows from
the same table. Our lift of the secondary element in direction `k` has its
pole at `180·k/N` degrees when `2k ≤ N` and at `180·(k−N)/N` otherwise
(the positive hemisphere of [Double groups](./double-groups.md), rotated by
`90°` under the twist); the reference's pole is the table entry. When the
two are antipodal the element is **swapped**: its reference lift is our
barred element. The `[barred.<group>]` tables list the remaining cases
explicitly (the polyhedral groups, whose poles the table does not cover).

### Polyhedral groups

The correspondence between canonical axis numbers and a reference's axis
names is a table, not a rule: `[overrides.<group>]` holds the complete map
of `T`, `Th`, `Td`, `O`, `Oh`, `I`, `Ih`.

## The Altmann & Herzig convention

Decided in wayfinder ticket
[Altmann & Herzig label-convention rules and overrides](https://github.com/hydeik/pointgroup-tables/issues/14)
from the transcriptions of
[Transcribe Altmann & Herzig element labels](https://github.com/hydeik/pointgroup-tables/issues/13).
The file `src/pointgroup_tables/data/conventions/labels/altmann-herzig.toml`
is generated by `tools/ah_convention.py generate` and never edited by hand.

- `exponent_style = "plusminus"`, `barred_format = "~{label}"`.
- `[secondary]` has the entries `M = 2, 4, …, 20`: `M ≤ 10` from the
  parameter tables of `D2h`, `D4h`, `D6h`, `D8h`, `D10h` (T 31.1 … T 39.1),
  `M = 12 … 20` from those of `D6d`, `D7h`, `D8d`, `D9h`, `D10d` (T 45.1,
  T 36.1, T 47.1, T 38.1, T 49.1), whose `C2'` and `sigma_v` / `sigma_d`
  rows are the `C2'` and `C2''` directions of `D_Mh`. The book's numbering
  is not a formula but its Rule 2 (subgroup consistency): `D9h` numbers the
  three axes of `D3h` first, `D8h` the two of `D4h` first, `D6d` its `C2'`
  as `C2'1 … C2'3, C2''1 … C2''3` of `D6h`. The book tabulates every axial
  family up to `n = 10`, so these entries cover its whole range; `n ≥ 11`
  is outside the book and gets canonical indices.
- `[overrides]`: the Cartesian names of `D2` (`C2z`, `C2x`, `C2y`), `D2h`
  (plus `sigma_z`, `sigma_x`, `sigma_y`) and `C2v` (`sigma_x`, `sigma_y`),
  and the complete maps of the seven polyhedral groups. `C2h`, `C2`, `Cs`,
  `Ci` need none.
- `[barred]`: `Td`, `O`, `Oh`, `I`, `Ih` (`T`, `Th` agree with the canonical
  section); the axial lists are derived.
- `[irreps]` as in [Irrep labels](./irrep-labels.md#reference-irrep-labels-in-a-label-convention),
  with the `[irreps.overrides]` the transcribed character tables decide
  (`REF-5`). The transcriptions cover the convention's whole declared range
  (every axial group with `n ≤ 10` and the seven polyhedral groups,
  [Reference irrep symbols of the `altmann-herzig` convention beyond the transcribed character tables](https://github.com/hydeik/pointgroup-tables/issues/36)),
  so every override is read off the book, none is extrapolated. What they
  hold: the `E_k` of `Dnd` with even `n` and of `S4n`, which the book
  numbers by the rotation `C_n` and the canonical rule by the principal
  element `S_2n` (`k' = k (n + 1) mod 2n`, folded by `E_k ≡ E_{2n−k}`:
  `D4d` `E1 ↔ E3`, `D6d` `E1 ↔ E5`, `D8d` `E1 ↔ E7`, `E3 ↔ E5`, `D10d`
  `E1 ↔ E9`, `E3 ↔ E7`, and `S8`, `S12`, `S16`, `S20` likewise); the
  members of the complex pairs of the cyclic-type groups (`Cn`, `Cnh`,
  `S2n`), where the book's `1E_k` is the canonical `2E_k` for every ordinary
  pair of odd `n` and for part of them at even `n` (`C8`, `C10` swap `E1`,
  `E2` and keep `E3`, `E4`), and of the projective pairs exactly those with
  `k ≡ 3 (mod 4)` in `E_{k/2}` except in `C4`, `C4h`, `S4`, which swap
  `1E1/2` instead; and the `F` spelling of the pairs of `T` and `Th`. The
  complex pairs follow no rule the specification states: the shipped file
  lists what the book prints, group by group. `D5h` needs no override:
  its printed character table (T 34.4) carries the symbols `E3/2` and
  `E7/2` on each other's rows, a recorded discrepancy that the book's own
  bases (T 34.6) and matrices (T 34.7) decide, corrected before the
  comparison. No `[classes.overrides]`.

**Element isomorphism.** The book's tables and the canonical labels are
matched by an isomorphism of abstract groups fixed by anchors (`E`, the
principal rotation, the first secondary element, `i`). For `C2v` and
`C6v` the anchor is the letter-preserving one (`sigma_y = sigma_v(1)`,
`sigma_v1 = sigma_v(1)`), not the geometric one (`sigma_x`, `sigma_d1`):
the book's `B1` of `C2v` is symmetric under `sigma_y` (T 50.4) and that of
`C6v` under `sigma_v` (Table 17.6), so only this choice lets rule 5 of the
irrep labels reproduce the book's symbols; it costs one extra swapped lift
in each group. For the polyhedral groups the anchors of the transcription
(`C4z+ = C4^1(1)`, `C3,1+ = C3^1(1)`, `C5,1+ = C5^1(1)`, second generator by
search) are kept: `T` and `Th` match the canonical senses and section
exactly, and for `O` / `I` no inner automorphism removes the sign
disagreements (best: 4 of 8 in `O`, 14 of 30 in `I`), which are the
book's axis orientations and pure data.

**Verification** (`REF-6`, `REF-4` of [Verification](./verification.md)):
`uv run python tools/ah_convention.py check` applies the rules and the
data to the canonical labels of every transcribed group and compares with
the book's labels element by element, and compares the derived swapped
list of every axial group with the list the transcription checker solves
from the book's tilde marks. Result at the time of the decision: all 33
axial files and both `C3v` settings reproduce with zero mismatches, the
derived swapped lists equal the solved ones in every group, and the
polyhedral tables reproduce by construction. The entries `M = 12 … 20`
have no transcribed group behind them and are checked only against the
printed parameter tables.

## Version 1 scope

The spec ships the canonical rules, the convention framework and the
complete Altmann & Herzig convention for elements. The book's `1E` / `2E`,
`B1` / `B2` and `T1` / `T2` assignments beyond `C2v` and `C6v` are checked
by `REF-5` once the character tables are transcribed
([Transcribe Altmann & Herzig character tables and sample axial irrep matrices](https://github.com/hydeik/pointgroup-tables/issues/19)).
Bradley & Cracknell and Bilbao conventions are later versions
([Command line](./cli.md#not-in-version-1)).
