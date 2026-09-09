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
| `barred` | `true` or `false`; double-group tables only |
| `base_index` | index of the point-group element this element covers; double-group tables only |

## Label conventions (reference labels)

One TOML file per convention, shipped with the package. The default is the
canonical labelling itself (no file).

```toml
name = "altmann-herzig"
source = "S. L. Altmann and P. Herzig, Point-Group Theory Tables, 2nd ed."
exponent_style = "plusminus"      # C3^1 -> C3+, C3^-1 -> C3-, C8^3 -> C8^3+
sigma_spelling = { sigma_h = "sigma_h", sigma_v = "sigma_v", sigma_d = "sigma_d" }

[rules]                           # applied to every group
swap_primes = false               # C2' <-> C2'' and sigma_v <-> sigma_d (even n)
axis_index_shift = 0              # rotate secondary-axis numbering

[overrides.D4h]                   # per-group explicit map, canonical -> reference
"sigma_v(1)" = "sigma_x"

barred_format = "~{label}"        # how the reference spells a barred element; "{label}~" for a suffix

[barred.D3h]                      # per-group: elements whose reference lift is our barred element
swapped = ["C2'(2)", "C2'(3)"]
```

The rule vocabulary of version 1 is exactly `exponent_style`,
`sigma_spelling`, `swap_primes` and `axis_index_shift`; anything a reference
does that these cannot express goes in `overrides`. A group without an
override table still receives rule-based reference labels. `barred_format`
and the `barred` tables concern the double group: they rename the bar and
swap the section for the listed elements, which changes the factor table and
the signs of the projective irreps as described in
[Double groups](./double-groups.md). Version 1 ships no `barred` table.

## Command line and metadata

- `--label-convention <name>` selects a convention (independent of
  `--basis-convention`, which selects the basis of irrep matrices).
- Output metadata records `label_convention = { name, source }` when one is
  selected.

## Version 1 scope

The spec ships the canonical rules and the convention framework. The
Altmann & Herzig convention starts with rules only and an empty override
table; filling it is the task
[Transcribe Altmann & Herzig element labels](https://github.com/hydeik/pointgroup-tables/issues/13)
followed by
[Altmann & Herzig label-convention rules and overrides](https://github.com/hydeik/pointgroup-tables/issues/14).
Bradley & Cracknell and Bilbao conventions are not yet specified.
