# Values

Decided in wayfinder ticket
[Serialization of exact cyclotomic values and their numeric approximations](https://github.com/hydeik/pointgroup-tables/issues/6).
Implements ADR-0002: every number in every table is exact, and floating
point is only an emitted approximation. The arithmetic is the tool's own
([ADR-0009](../adr/0009-no-runtime-dependencies.md)): no SymPy, no mpmath.

## The field of a document

Every value in a document is an element of one cyclotomic field `Q(ζ_M)`,
`ζ_M = exp(2πi/M)`, declared once in the document metadata:

```json
"field": { "kind": "cyclotomic", "conductor": 12, "degree": 4,
           "generator": "exp(2*pi*i/12)", "basis": "power" }
```

`M` (the **conductor**) is fixed per group and basis convention: it is the
smallest `M` such that every entry of every irrep matrix of that group, in
that basis convention, lies in `Q(ζ_M)`. It does not depend on which table
subcommand produced the document, so a given value always has the same
encoding whether it appears in `characters`, `matrices` or `all`. Values
computed internally in a smaller field `Q(ζ_m)` are embedded at the output
boundary through `ζ_m = ζ_M^{M/m}`. Minimality is a property of the
**unfiltered** document (`--irreps all`, both kinds of irrep): a filtered
document may fit in a smaller field and still declares the same `M`, so
that the field never depends on a filter (`VL-3` of
[Verification](./verification.md)).

### Chain and coupling documents

Version 2 adds two cases, both fixed by the group and the options alone
and declared in the same `field` record without new keys:

- A document produced with `--chain` declares the **chain conductor**
  `M_chain` of its group and chain class
  ([Group chains](./group-chains.md#the-conductor-of-a-chain-document)):
  the smallest multiple of the chain field `lcm(M_G, M_{H_i})` in which
  every chain-adapted matrix of the group can be written, which is the
  version-1 `M_G` for every chain of the axial families (so the
  default-chain document of an axial group has the version-1 field and
  values) and a per-chain value for the polyhedral groups (`8` for
  `O ⊃ C4`, `12` for `T ⊃ C3`, `24` for `T ⊃ D2 ⊃ C2`, `20` for
  `I ⊃ D5 ⊃ C5`, `60` for `I ⊃ T ⊃ C3`, …; `list chains` prints it). The
  values of the chain-adapted matrices and of the `subgroup_chain` table
  are elements of this field through the canonical norm factor of
  [Exact computation](./group-chains.md#exact-computation): no square
  root outside the field is ever taken.
- A `cg` document declares `lcm(M_basis, 8)` for the axial families and
  `lcm(M_basis, 24)` for the polyhedral groups, `M_basis` being the
  conductor a `matrices` document with the same options declares
  (`M_G`, the convention's, or `M_chain`)
  ([Clebsch–Gordan coefficients](./clebsch-gordan.md#the-field-of-a-cg-document));
  `cg` is the one subcommand whose field differs from that of the other
  tables of the same group and options.

`degree` is `φ(M)`; the basis is the power basis `1, ζ, ζ², ..., ζ^{φ(M)−1}`.

## Exact form

An exact value is the dense coefficient vector of the element over the power
basis, each coefficient a rational number written as an integer pair
`[numerator, denominator]`:

```
exact value := [ [p_0, q_0], [p_1, q_1], ..., [p_{φ(M)−1}, q_{φ(M)−1}] ]
```

Rules:

- exactly `φ(M)` entries, in ascending power of `ζ`;
- `q_i > 0` and `gcd(p_i, q_i) = 1`; zero is `[0, 1]`;
- the vector is the canonical reduced representative (degree `< φ(M)`), so
  two values are equal exactly when their vectors are equal.

Example, `M = 3` (`φ = 2`, basis `1, ζ`): `ζ` is `[[0,1],[1,1]]`, `ζ⁻¹ = ζ² =
−1 − ζ` is `[[-1,1],[-1,1]]`, `1/2` is `[[1,2],[0,1]]`.

No strings are involved: any language with integers can read an exact value,
and a Python consumer rebuilds it as `sum(Fraction(p, q) * zeta**i)`.

## Numeric form

A numeric value is always a pair of doubles `[re, im]`, including for real
values (`[x, 0.0]`). It is obtained by evaluating the exact value with at
least 30 significant digits and rounding once to double precision, so it is
the correctly rounded approximation. The evaluation uses the standard
library's `decimal` at 40 significant digits: `cos(2πj/M)` and
`sin(2πj/M)` from their Taylor series (with `π` to the same precision), the
rational coefficients exactly, the sum in `decimal`, and one final
conversion of each part to `float` (round half to even). It is printed with
the shortest representation that round-trips (Python `repr(float)`, at
most 17 significant digits); negative zero is printed as `0.0`. The test
suite compares the values with an independent high-precision evaluation
(mpmath as a development dependency) on every group of the `seed` set.

## Pretty form (optional)

A pretty value is a string in the value grammar of
[Basis conventions](./basis-conventions.md#value-grammar), such as
`"-1 - zeta"`, `"zeta^7"`, `"(1 + sqrt5)/4"` or `"-1/2*i*sqrt3"`, with
`zeta` meaning `ζ_M`; the tool writes the shortest of a rational,
`r*zeta^k`, `r*sqrtN*zeta^k` and `re + i*im` with `re`, `im` rational
combinations of `1`, `sqrt2`, `sqrt3`, `sqrt5` and their products, as the
reference transcriptions do. When none of these spellings exists (the
chain-adapted matrices of `I ⊃ T ⊃ C3` and their norm factors are elements
of degree eight of `Q(ζ_60)⁺` whose real parts involve `cos(2π/15)`), the
pretty form is the polynomial in `zeta` over the power basis,
`p_0 + p_1*zeta + p_2*zeta^2 + …`, with the rational coefficients of the
exact form; no further symbol is introduced. The grammar parses exactly,
but the spelling of a value is not unique, so the pretty form is for human
readers only and must not be used for comparison or reconstruction. It is
emitted only on request.

## Layout in tables

Every table that holds values carries parallel arrays of identical shape,
one per requested form:

```json
{ "exact":   [[[[0,1],[1,1]], [[0,1],[0,1]]],
              [[[0,1],[0,1]], [[-1,1],[-1,1]]]],
  "numeric": [[[-0.5, 0.8660254037844386], [0.0, 0.0]],
              [[0.0, 0.0], [-0.5, -0.8660254037844386]]] }
```

(the matrix `diag(ζ, ζ⁻¹)` of `C3` in the `E` irrep of `C3v`, `M = 3`).

- Matrices are nested arrays, rows then columns, without a separate `shape`.
  The numeric form of a `d × d` matrix is therefore a `[d][d][2]` tensor and
  the exact form a `[d][d][φ(M)][2]` tensor of integers.
- A character table is laid out the same way with irreps as rows and
  conjugacy classes as columns.
- The forms present are chosen by `--values exact|numeric|both|pretty`, a
  comma-separated set; the default is `both`. A form that is not requested is
  omitted, never `null`.

## YAML and TOML

YAML and TOML documents are direct dumps of the same document model. Nested
arrays are valid in both; TOML has no `null`, and no `null` ever appears in
the model (absent keys instead). No format-specific layout exists.
