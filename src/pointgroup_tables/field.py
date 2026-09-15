"""The document field ``Q(zeta_M)`` and its exact arithmetic
(docs/spec/values.md, docs/spec/basis-conventions.md, ADR-0009).

An element of ``Q(zeta_M)`` is a tuple of ``phi(M)`` rationals, its
coefficient vector over the power basis ``1, zeta, ..., zeta^(phi(M) - 1)``
reduced modulo the cyclotomic polynomial. Two elements are equal exactly when
their tuples are equal. The arithmetic is the prototypes' promoted to
production: no SymPy, no mpmath.
"""

from __future__ import annotations

import math
from fractions import Fraction
from functools import cache
from typing import Any

from .rotation import RotationGroup

Value = tuple[Fraction, ...]
Matrix = list[list[Value]]

POLYHEDRAL_CONDUCTOR = {"T": 12, "O": 8, "I": 20}


def euler_phi(m: int) -> int:
    return sum(1 for k in range(1, m + 1) if math.gcd(k, m) == 1)


def canonical_conductor(R: RotationGroup) -> int:
    """The canonical conductor, fixed by ``R`` alone."""
    if R.axis_order is None:
        return POLYHEDRAL_CONDUCTOR[R.name]
    N = R.axis_order
    if R.name.startswith("D"):
        return math.lcm(2 * N, 4)
    return N if N % 2 else 2 * N


def field_metadata(conductor: int) -> dict[str, Any]:
    return {
        "kind": "cyclotomic",
        "conductor": conductor,
        "degree": euler_phi(conductor),
        "generator": f"exp(2*pi*i/{conductor})",
        "basis": "power",
    }


# ---- cyclotomic polynomials ----------------------------------------------------


def _divide_exactly(numerator: list[int], denominator: list[int]) -> list[int]:
    """Exact division of integer polynomials (coefficients low to high)."""
    rest = list(numerator)
    quotient = [0] * (len(rest) - len(denominator) + 1)
    for k in range(len(rest) - len(denominator), -1, -1):
        c = rest[k + len(denominator) - 1] // denominator[-1]
        quotient[k] = c
        for j, d in enumerate(denominator):
            rest[k + j] -= c * d
    if any(rest):
        raise ArithmeticError("the polynomial division is not exact")
    return quotient


@cache
def cyclotomic_polynomial(m: int) -> tuple[int, ...]:
    """``Phi_m`` as integer coefficients, low to high, monic."""
    poly = [-1] + [0] * (m - 1) + [1]  # x^m - 1
    for d in range(1, m):
        if m % d == 0:
            poly = _divide_exactly(poly, list(cyclotomic_polynomial(d)))
    return tuple(poly)


# ---- the field ------------------------------------------------------------------


class CyclotomicField:
    """``Q(zeta_M)`` with ``zeta_M = exp(2 pi i / M)``; values are reduced tuples."""

    def __init__(self, conductor: int) -> None:
        if conductor < 1:
            raise ValueError("the conductor must be positive")
        self.conductor = conductor
        self.polynomial = cyclotomic_polynomial(conductor)
        self.degree = len(self.polynomial) - 1
        # zeta^k for k in [degree, 2 degree - 2] as reduced vectors: one table
        # lookup per high power turns the reduction into a linear combination
        self._high: list[Value] = []
        for k in range(self.degree, 2 * self.degree - 1):
            self._high.append(self._reduce_power(k))
        self._zeta: list[Value] = [self._reduce_power(k) for k in range(conductor)]

    def __repr__(self) -> str:
        return f"CyclotomicField({self.conductor})"

    def __eq__(self, other: object) -> bool:
        return isinstance(other, CyclotomicField) and other.conductor == self.conductor

    def __hash__(self) -> int:
        return hash(("CyclotomicField", self.conductor))

    def _reduce_power(self, k: int) -> Value:
        """``zeta^k`` reduced by long division (used only to build the tables)."""
        n = self.degree
        coefficients = [Fraction(0)] * (k + 1)
        coefficients[k] = Fraction(1)
        for j in range(k, n - 1, -1):
            q = coefficients[j]
            if q:
                for t, p in enumerate(self.polynomial):
                    coefficients[j - n + t] -= q * p
        coefficients = coefficients[:n] + [Fraction(0)] * (n - len(coefficients))
        return tuple(coefficients[:n])

    # ---- construction ---------------------------------------------------------

    def reduce(self, coefficients: list[Fraction]) -> Value:
        """Reduce a coefficient list of any length (low to high) to a value."""
        n = self.degree
        out = list(coefficients[:n]) + [Fraction(0)] * max(0, n - len(coefficients))
        for k in range(n, len(coefficients)):
            a = coefficients[k]
            if a:
                if k - n >= len(self._high):
                    raise ValueError("a power beyond the reduction table")
                for j, h in enumerate(self._high[k - n]):
                    if h:
                        out[j] += a * h
        return tuple(out)

    def rational(self, r: Fraction | int) -> Value:
        return (Fraction(r),) + (Fraction(0),) * (self.degree - 1)

    @property
    def zero(self) -> Value:
        return self.rational(0)

    @property
    def one(self) -> Value:
        return self.rational(1)

    def zeta(self, k: int = 1) -> Value:
        """``zeta_M^k``."""
        return self._zeta[k % self.conductor]

    def root_of_unity(self, m: int, k: int = 1) -> Value:
        """``zeta_m^k`` as an element of this field.

        ``zeta_m`` lies in ``Q(zeta_M)`` when ``m | M``, and also when ``m = 2 m'``
        with ``m'`` odd and ``m' | M`` (``zeta_2m' = -zeta_m'^((m' + 1) / 2)``).
        """
        M = self.conductor
        if M % m == 0:
            return self.zeta((M // m) * k)
        if m % 2 == 0 and (m // 2) % 2 == 1 and M % (m // 2) == 0:
            half = m // 2
            step = M // half
            # zeta_2h = -zeta_h^((h + 1) / 2); its k-th power is (-1)^k zeta_h^(k (h + 1) / 2)
            value = self.zeta(step * (k * ((half + 1) // 2)))
            return self.neg(value) if k % 2 else value
        raise ValueError(f"zeta_{m} does not lie in Q(zeta_{M})")

    def has_root_of_unity(self, m: int) -> bool:
        M = self.conductor
        return M % m == 0 or (m % 2 == 0 and (m // 2) % 2 == 1 and M % (m // 2) == 0)

    def constant(self, name: str) -> Value:
        """The named constants of the value grammar (docs/spec/basis-conventions.md)."""
        M = self.conductor
        if name == "i":
            if M % 4:
                raise ValueError(f"'i' needs 4 | M, the conductor is {M}")
            return self.zeta(M // 4)
        if name == "sqrt2":
            if M % 8:
                raise ValueError(f"'sqrt2' needs 8 | M, the conductor is {M}")
            return self.add(self.zeta(M // 8), self.zeta(-(M // 8)))
        if name == "sqrt3":
            if M % 12:
                raise ValueError(f"'sqrt3' needs 12 | M, the conductor is {M}")
            return self.add(self.zeta(M // 12), self.zeta(-(M // 12)))
        if name == "sqrt5":
            if M % 5:
                raise ValueError(f"'sqrt5' needs 5 | M, the conductor is {M}")
            z = M // 5
            return self.sub(
                self.add(self.zeta(z), self.zeta(4 * z)),
                self.add(self.zeta(2 * z), self.zeta(3 * z)),
            )
        if name == "phi":
            return self.scale(
                self.add(self.one, self.constant("sqrt5")), Fraction(1, 2)
            )
        raise ValueError(f"unknown constant '{name}'")

    def has_constant(self, name: str) -> bool:
        M = self.conductor
        needs = {"i": 4, "sqrt2": 8, "sqrt3": 12, "sqrt5": 5, "phi": 5}
        return M % needs[name] == 0

    # ---- arithmetic -----------------------------------------------------------

    @staticmethod
    def add(x: Value, y: Value) -> Value:
        return tuple(a + b for a, b in zip(x, y, strict=True))

    @staticmethod
    def sub(x: Value, y: Value) -> Value:
        return tuple(a - b for a, b in zip(x, y, strict=True))

    @staticmethod
    def neg(x: Value) -> Value:
        return tuple(-a for a in x)

    @staticmethod
    def scale(x: Value, r: Fraction | int) -> Value:
        return tuple(a * r for a in x)

    def mul(self, x: Value, y: Value) -> Value:
        n = self.degree
        product = [Fraction(0)] * (2 * n - 1)
        for i, a in enumerate(x):
            if a:
                for j, b in enumerate(y):
                    if b:
                        product[i + j] += a * b
        return self.reduce(product)

    def power(self, x: Value, k: int) -> Value:
        if k < 0:
            return self.power(self.inverse(x), -k)
        result = self.one
        base = x
        while k:
            if k & 1:
                result = self.mul(result, base)
            k >>= 1
            if k:
                base = self.mul(base, base)
        return result

    def conjugate(self, x: Value) -> Value:
        """Complex conjugation ``zeta -> zeta^-1``."""
        out = [Fraction(0)] * self.degree
        for j, a in enumerate(x):
            if a:
                for t, z in enumerate(self.zeta(-j)):
                    if z:
                        out[t] += a * z
        return tuple(out)

    def inverse(self, x: Value) -> Value:
        """``1 / x`` by solving the linear system ``x * y = 1`` over the rationals."""
        if self.is_zero(x):
            raise ZeroDivisionError("division by zero in the field")
        n = self.degree
        columns = [self.mul(x, self.zeta(j)) for j in range(n)]
        rows = [
            [columns[j][r] for j in range(n)] + [Fraction(int(r == 0))]
            for r in range(n)
        ]
        for col in range(n):
            pivot = next(r for r in range(col, n) if rows[r][col] != 0)
            rows[col], rows[pivot] = rows[pivot], rows[col]
            p = rows[col][col]
            rows[col] = [v / p for v in rows[col]]
            for r in range(n):
                if r != col and rows[r][col] != 0:
                    f = rows[r][col]
                    rows[r] = [
                        a - f * c for a, c in zip(rows[r], rows[col], strict=True)
                    ]
        return tuple(rows[r][n] for r in range(n))

    def divide(self, x: Value, y: Value) -> Value:
        return self.mul(x, self.inverse(y))

    def total(self, values: Any) -> Value:
        out = list(self.zero)
        for value in values:
            for j, a in enumerate(value):
                if a:
                    out[j] += a
        return tuple(out)

    # ---- predicates -----------------------------------------------------------

    @staticmethod
    def is_zero(x: Value) -> bool:
        return not any(x)

    @staticmethod
    def is_rational(x: Value) -> bool:
        return not any(x[1:])

    def is_real(self, x: Value) -> bool:
        return self.conjugate(x) == x

    def is_imaginary(self, x: Value) -> bool:
        """Purely imaginary (or zero)."""
        return self.conjugate(x) == self.neg(x)

    @staticmethod
    def as_rational(x: Value) -> Fraction:
        if any(x[1:]):
            raise ValueError("the value is not rational")
        return x[0]

    def norm_squared(self, x: Value) -> Value:
        return self.mul(x, self.conjugate(x))

    # ---- embeddings ------------------------------------------------------------

    def lift(self, x: Value, target: CyclotomicField) -> Value:
        """Embed a value into ``Q(zeta_M')`` with ``M | M'`` through ``zeta_M = zeta_M'^(M'/M)``."""
        if target.conductor % self.conductor:
            raise ValueError(
                f"Q(zeta_{self.conductor}) does not embed by conductor in Q(zeta_{target.conductor})"
            )
        if target.conductor == self.conductor:
            return x
        step = target.conductor // self.conductor
        return target.total(
            target.scale(target.zeta(j * step), a) for j, a in enumerate(x) if a
        )

    def descend(self, x: Value, target: CyclotomicField) -> Value | None:
        """The same value in the subfield ``Q(zeta_M')`` with ``M' | M``, or ``None``
        when it does not lie there."""
        if self.conductor % target.conductor:
            raise ValueError("the target conductor does not divide this one")
        if target.conductor == self.conductor:
            return x
        n = target.degree
        columns = [target.lift(target.zeta(j), self) for j in range(n)]
        solution = _solve_rational(columns, x)
        return None if solution is None else tuple(solution)

    # ---- output forms ----------------------------------------------------------

    @staticmethod
    def exact(x: Value) -> list[list[int]]:
        """The exact form: ``phi(M)`` integer pairs ``[numerator, denominator]``."""
        return [[a.numerator, a.denominator] for a in x]

    def from_exact(self, pairs: list[list[int]]) -> Value:
        if len(pairs) != self.degree:
            raise ValueError(f"an exact value needs {self.degree} coefficients")
        return tuple(Fraction(p, q) for p, q in pairs)


def _solve_rational(columns: list[Value], target: Value) -> list[Fraction] | None:
    """Rational coordinates of ``target`` over the given vectors, or ``None``."""
    rows = len(target)
    cols = len(columns)
    matrix = [[columns[j][r] for j in range(cols)] + [target[r]] for r in range(rows)]
    pivots: list[int] = []
    row = 0
    for col in range(cols):
        pivot = next((r for r in range(row, rows) if matrix[r][col] != 0), None)
        if pivot is None:
            continue
        matrix[row], matrix[pivot] = matrix[pivot], matrix[row]
        p = matrix[row][col]
        matrix[row] = [v / p for v in matrix[row]]
        for r in range(rows):
            if r != row and matrix[r][col] != 0:
                f = matrix[r][col]
                matrix[r] = [
                    a - f * c for a, c in zip(matrix[r], matrix[row], strict=True)
                ]
        pivots.append(col)
        row += 1
    if any(matrix[r][cols] != 0 for r in range(row, rows)):
        return None
    solution = [Fraction(0)] * cols
    for r, col in enumerate(pivots):
        solution[col] = matrix[r][cols]
    return solution


@cache
def cyclotomic_field(conductor: int) -> CyclotomicField:
    """The field of a conductor, built once."""
    return CyclotomicField(conductor)


# ---- matrices over the field -------------------------------------------------------


def matrix_identity(F: CyclotomicField, d: int) -> Matrix:
    return [[F.one if i == j else F.zero for j in range(d)] for i in range(d)]


def matrix_mul(F: CyclotomicField, A: Matrix, B: Matrix) -> Matrix:
    n, m, p = len(A), len(B), len(B[0])
    return [
        [
            F.total(F.mul(A[i][k], B[k][j]) for k in range(m) if not F.is_zero(A[i][k]))
            for j in range(p)
        ]
        for i in range(n)
    ]


def matrix_power(F: CyclotomicField, A: Matrix, k: int) -> Matrix:
    result = matrix_identity(F, len(A))
    for _ in range(k):
        result = matrix_mul(F, result, A)
    return result


def matrix_neg(F: CyclotomicField, A: Matrix) -> Matrix:
    return [[F.neg(a) for a in row] for row in A]


def matrix_trace(F: CyclotomicField, A: Matrix) -> Value:
    return F.total(A[i][i] for i in range(len(A)))
