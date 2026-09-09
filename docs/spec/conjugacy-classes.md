# Conjugacy classes

Decided in wayfinder ticket
[Conjugacy classes: labels, representatives and projective (double-group) classes](https://github.com/hydeik/pointgroup-tables/issues/9).
Builds on [Abstract groups](./abstract-groups.md) (canonical element order),
[Label conventions](./label-conventions.md) (canonical element labels),
[Double groups](./double-groups.md) (`G†`, the canonical section `s`, the
element order of `G†`) and [Irrep labels](./irrep-labels.md) (the rows of
the character table). Classes of any group can be printed with
`prototypes/abstract_group_prototype.py labels <symbols>`.

## Classes of `G`

Conjugacy classes are computed exactly from the Cayley table. They are
**ordered by their first member** in the canonical element order, and the
**representative** of a class is that first member. Because the canonical
order lists proper elements before improper ones and keeps the classes of
`R` contiguous, the classes of `D6h` come out as

```
E  2C6  2C3  C2  3C2'  3C2''  i  2S3  2S6  sigma_h  3sigma_d  3sigma_v
```

which is the order of Altmann & Herzig. The order is not always theirs
(`D3h`: `E 2C3 3C2' 2S3 sigma_h 3sigma_v`, the book puts `sigma_h` before
`2S3`); the class index is canonical and **a label convention never
reorders classes**. A consumer that wants a reference order permutes
columns by index.

The representative of every class used by the irrep rules is the element
those rules name: `C4^1(1)` is the first member of the `6C4` class of `Oh`,
`C5^1(1)` the first member of `12C5` in `I`, the principal and secondary
elements `(1, 0)` and `(0, 1)` the first members of their classes.

## Class labels

A **class label** is the size of the class followed by a reduced form of
the representative's canonical label. ASCII, no whitespace:

```
class_label := block ("+" block)?        # two blocks only for a merged double-group class
block       := size? "~"? rep
size        := integer ≥ 2               # omitted when the class has one element
rep         := "E" | "i" | "C2" | "C2'" | "C2''"
             | "C" m ("^" j)? | "S" m ("^" j)?
             | "sigma_h" | "sigma_v" | "sigma_d" | "sigma"
```

`rep` is the canonical label of the representative with

1. the axis index `(k)` removed (`C2'(1)` → `C2'`, `sigma_d(2)` →
   `sigma_d`, `C3^1(4)` → `C3`);
2. the exponent `^1` removed always (`C6^1` → `C6`);
3. the sign of the exponent removed when the class contains the inverse of
   the representative (`{C3^1, C3^-1}` → `C3`, `{C5^2, C5^-2}` → `C5^2`,
   `{S10^3, S10^-3}` → `S10^3`); kept otherwise (`Th`: `4C3` and `4C3^-1`,
   `4S6^-1` and `4S6`; `C6`: `C6 C3 C2 C3^-1 C6^-1`).

Regular expression accepting exactly these labels:

```
^(\d*~?(E|i|C2''?|[CS]\d+(\^-?\d+)?|sigma(_[hvd])?))(\+\d*~(E|i|C2''?|[CS]\d+(\^-?\d+)?|sigma(_[hvd])?))?$
```

Class labels are unique within a group: two classes with the same reduced
representative would have the same kind, order, exponent modulus and
prime, and the canonical element rules give such elements the same axis
class only when they are conjugate. Examples: `Td`: `E 3C2 8C3 6sigma_d
6S4`; `Oh`: `E 3C2 6C2' 8C3 6C4 i 3sigma_h 6sigma_d 8S6 6S4`; `Ih`: `E
15C2 20C3 12C5 12C5^2 i 15sigma 20S6 12S10^3 12S10`; `Th`: `E 3C2 4C3
4C3^-1 i 3sigma_h 4S6^-1 4S6`; `D2d`: `E C2 2C2' 2S4 2sigma_d`.

## Classes of `G†`

Let `K` be a class of `G` with representative `g`. Its preimage in `G†` is
either two classes, the one containing `s(g)` and the one containing
`Ē·s(g)` (`K` **splits**), or a single class of size `2|K|` (`K` is
**merged**). `K` splits exactly when it is **ω-regular**, that is when
projective characters can be non-zero on it; the number of regular classes
equals the number of projective irreps.

Under the canonical section a split class is not, in general, the set
`s(K)`: in `D3†` the class of `s(C2'(1))` is `{C2'(1), ~C2'(2), ~C2'(3)}`
and the class of `Ē·s(C2'(1))` is `{~C2'(1), C2'(2), C2'(3)}`, so the
character of `1E3/2` is `+i` on `C2'(1)` and `−i` on `C2'(2)`. Projective
characters are therefore class functions of `G†`, not of `G`, and every
table of projective characters is indexed by the classes of `G†`. (Altmann
& Herzig bar `C2'(2)` and `C2'(3)` instead, which makes their `3C2'` a
single class; that is the `barred` swap of the label convention described
in [Double groups](./double-groups.md), and it changes the class membership
accordingly.)

The classes of `G†` are ordered **by base class**: for every class `K` of
`G` in canonical class order, first the class containing `s(g)`, then, if
`K` splits, the class containing `Ē·s(g)`. Their representatives are
`s(g)` and `Ē·s(g)` (the element labelled `~` plus the label of `g`). The
labels follow the grammar above:

- a split class: `size` in `G†` plus the reduced representative, the barred
  one with `~` (`E` and `~E`, `2C3` and `2~C3`, `3C2'` and `3~C2'`);
- a merged class: the two blocks joined by `+`, each with the size of its
  half (`C2+~C2`, `2C2'+2~C2'`, `3sigma_v+3~sigma_v`).

Examples (class counts as in [Double groups](./double-groups.md)):

| Group | Classes of `G†` |
| ----- | --------------- |
| `C6†` | `E ~E C6 ~C6 C3 ~C3 C2 ~C2 C3^-1 ~C3^-1 C6^-1 ~C6^-1` |
| `D3†` | `E ~E 2C3 2~C3 3C2' 3~C2'` |
| `D4†` | `E ~E 2C4 2~C4 C2+~C2 2C2'+2~C2' 2C2''+2~C2''` |
| `D3h†` | `E ~E 2C3 2~C3 3C2'+3~C2' 2S3 2~S3 sigma_h+~sigma_h 3sigma_v+3~sigma_v` |
| `T†` | `E ~E 3C2+3~C2 4C3 4~C3 4C3^-1 4~C3^-1` |
| `O†` | `E ~E 3C2+3~C2 6C2'+6~C2' 8C3 8~C3 6C4 6~C4` |

## Class records

Classes are part of the **`elements` table** (there is no `classes`
subcommand, see [Command line](./cli.md)); with `--double-group` they are
the classes of `G†`. Every element record gains the field `class`, the
index of its class. Each class is emitted as:

| Field | Values |
| ----- | ------ |
| `index` | position in the canonical class order |
| `label` | canonical class label |
| `reference_label` | present only when a label convention is selected |
| `size` | number of elements |
| `representative` | element index of the representative |
| `members` | element indices in canonical order (the representative first) |
| `inverse` | index of the class of the inverses of the members (itself for an ambivalent class) |
| `regular` | `G` only: `true` when the class splits in `G†` |
| `base_class` | `G†` only: index of the class of `G` the members cover |
| `split` | `G†` only: `true` for a split class, `false` for a merged one |
| `barred` | `G†` only, split classes only: `true` for the class containing `Ē·s(g)` |

`inverse` links the classes whose characters are complex conjugates of
each other; `regular` counts the projective irreps.

## Columns of the character table

The `characters` table is indexed by class index and declares which group
its columns belong to:

| `--irreps` | Columns | `class_group` |
| ---------- | ------- | ------------- |
| `ordinary` | classes of `G` | `"G"` |
| `projective`, `all` | classes of `G†` (ordered by base class as above) | `"G_double"` |

With `all` the table has one row per ordinary and projective irrep and one
column per class of `G†`; an ordinary irrep repeats its value on the barred
class of a split pair, a projective irrep negates it. This is the layout of
the double-group tables of Altmann & Herzig. Whether the class list is
repeated inside the `characters` table or only referenced by index is
decided by the output-schema chapter.

## Reference class labels

When a label convention is selected, the reference class label is derived
from the **reference label of the representative** by the same reduction
(axis index removed, exponent sign removed when the class contains the
inverse), so with `exponent_style = "plusminus"` the `Th` classes read
`4C3+`, `4C3-`, `4S6-`, `4S6+` and the `D6h` class `{C3+, C3-}` reads
`2C3`. A convention may add `[classes.overrides.<group>]` mapping canonical
class labels to the reference's spelling where the reference is irregular;
the convention cannot change the class order or membership except through
its `barred` table.

## Verification

An implementation must check, for every group it constructs:

- the classes partition the elements, class sizes divide the group order,
  and `size × |centraliser|` equals the order for the representative;
- the class labels match the grammar and are pairwise distinct;
- the class counts of `G` and `G†` agree with the table in
  [Double groups](./double-groups.md), and the number of `regular` classes
  equals the number of projective irreps;
- `inverse` is an involution on class indices; `base_class`, `split` and
  `barred` are consistent with the `base_index` and `barred` fields of the
  element records;
- for the polyhedral groups the class labels reproduce the class structures
  listed in [Label conventions](./label-conventions.md).

Orthogonality of the character table over these classes belongs to the
verification chapter.
