# Values

Decided in wayfinder ticket
[Serialization of exact cyclotomic values and their numeric approximations](https://github.com/hydeik/pointgroup-tables/issues/6).
Implements ADR-0002: every number in every table is exact, and floating
point is only an emitted approximation.

## The field of a document

Every value in a document is an element of one cyclotomic field `Q(ζ_M)`,
`ζ_M = exp(2πi/M)`, declared once in the document metadata:

```json
"field": { "kind": "cyclotomic", "conductor": 12, "degree": 4,
           "generator": "exp(2*pi*i/12)", "basis": "power" }
```

`M` (the **conductor**) is fixed per group and basis convention: it is the
smallest `M` such that every entry of every irrep matrix of that group, in
that basis convention, lies in `Q(ζ_M)`. It does not depend on which
subcommand produced the document, so a given value always has the same
encoding whether it appears in `characters`, `matrices` or `all`. Values
computed internally in a smaller field `Q(ζ_m)` are embedded at the output
boundary through `ζ_m = ζ_M^{M/m}`.

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
least 30 significant digits (mpmath) and rounding once to double precision,
so it is the correctly rounded approximation. It is printed with the shortest
representation that round-trips (Python `repr(float)`, at most 17 significant
digits); negative zero is printed as `0.0`.

## Pretty form (optional)

A pretty value is a simplified SymPy-style string such as `"-1 - zeta"` or
`"(1 + sqrt(5))/4"`, with `zeta` meaning `ζ_M`. Simplified forms are not
unique, so the pretty form is for human readers only and must not be used for
comparison or reconstruction. It is emitted only on request.

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
