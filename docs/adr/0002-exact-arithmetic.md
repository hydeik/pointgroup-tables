---
status: accepted
---

# All computation is exact; floating-point values are emitted, never computed with

Characters and irrep matrix entries live in cyclotomic fields, so every
internal computation (Cayley tables, factor systems, characters, irrep
matrices, orthogonality checks) is done exactly with SymPy. Double-precision
complex values appear only in the output, as an approximation printed next to
each exact value for consumers that cannot evaluate symbolic expressions.

## Considered Options

- **NumPy floating point throughout**: simple and fast, but equality of
  elements and characters becomes tolerance-based, and outputs cannot be
  verified against reference tables symbolically. Rejected.
- **Exact throughout with SymPy** (chosen): slower and adds a dependency, but
  every table is reproducible and provable, and the numeric output is derived
  from it in one place.

## Consequences

- Performance for large axis order `n` is a known cost, to be addressed only
  if it becomes a problem.
- The output format must carry both an exact and a numeric form of every
  value.
