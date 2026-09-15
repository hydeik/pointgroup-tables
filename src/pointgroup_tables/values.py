"""The three forms of a value (docs/spec/values.md): the exact coefficient
vector, the correctly rounded numeric pair and the optional pretty string in
the value grammar of docs/spec/basis-conventions.md.

The numeric form is evaluated with the standard library's ``decimal`` at
``PRECISION`` significant digits (``cos`` and ``sin`` from their Taylor
series, ``pi`` to the same precision) and rounded once to double. The pretty
form has one printer and one parser; the parser also reads the convention
data and the reference transcriptions.
"""

from __future__ import annotations

import re
from decimal import ROUND_HALF_EVEN, Context, Decimal
from fractions import Fraction
from functools import cache
from itertools import combinations

from .field import CyclotomicField, Value, _solve_rational

PRECISION = 40
# pi to 60 significant digits; the numeric evaluation rounds it to PRECISION
PI_DIGITS = "3.14159265358979323846264338327950288419716939937510582097494459"

# ---- numeric form ------------------------------------------------------------------


def _context() -> Context:
    return Context(prec=PRECISION, rounding=ROUND_HALF_EVEN)


@cache
def _pi() -> Decimal:
    return _context().create_decimal(PI_DIGITS)


def _cos_sin(theta: Decimal, ctx: Context) -> tuple[Decimal, Decimal]:
    """``cos`` and ``sin`` of ``theta`` (``|theta| <= pi``) from the Taylor series."""
    cos_total, sin_total = Decimal(0), Decimal(0)
    term = Decimal(1)  # theta^k / k!
    k = 0
    epsilon = ctx.create_decimal(10) ** (-(PRECISION + 5))
    while True:
        if k % 2 == 0:
            cos_total = ctx.add(cos_total, term if k % 4 == 0 else -term)
        else:
            sin_total = ctx.add(sin_total, term if k % 4 == 1 else -term)
        k += 1
        term = ctx.divide(ctx.multiply(term, theta), Decimal(k))
        if abs(term) < epsilon:
            break
    return cos_total, sin_total


@cache
def _root_of_unity(j: int, M: int) -> tuple[Decimal, Decimal]:
    """``cos`` and ``sin`` of ``2 pi j / M`` with the angle reduced to ``(-pi, pi]``."""
    ctx = _context()
    j %= M
    if 2 * j > M:
        j -= M
    theta = ctx.divide(
        ctx.multiply(ctx.multiply(Decimal(2), _pi()), Decimal(j)), Decimal(M)
    )
    return _cos_sin(theta, ctx)


def _to_double(x: Decimal) -> float:
    value = float(x)  # round half to even
    return 0.0 if value == 0 else value


def numeric(F: CyclotomicField, x: Value) -> list[float]:
    """The numeric form ``[re, im]``: the correctly rounded doubles of the value.

    A rational value is rounded directly; a value that is exactly real (or
    exactly imaginary) has its other part set to ``0.0`` exactly, so that
    cancellation in the decimal sum never leaves a spurious residue.
    """
    if F.is_rational(x):
        return [
            _to_double(Decimal(x[0].numerator) / Decimal(x[0].denominator))
            if x[0]
            else 0.0,
            0.0,
        ]
    ctx = _context()
    re_total, im_total = Decimal(0), Decimal(0)
    M = F.conductor
    for j, a in enumerate(x):
        if a:
            coefficient = ctx.divide(Decimal(a.numerator), Decimal(a.denominator))
            cos_j, sin_j = _root_of_unity(j, M)
            re_total = ctx.add(re_total, ctx.multiply(coefficient, cos_j))
            im_total = ctx.add(im_total, ctx.multiply(coefficient, sin_j))
    re_value = 0.0 if F.is_imaginary(x) else _to_double(re_total)
    im_value = 0.0 if F.is_real(x) else _to_double(im_total)
    return [re_value, im_value]


# ---- the value grammar: printer ----------------------------------------------------

_ROOTS = (("sqrt2", 2, 8), ("sqrt3", 3, 12), ("sqrt5", 5, 5))


def format_rational(r: Fraction) -> str:
    return str(r.numerator) if r.denominator == 1 else f"{r.numerator}/{r.denominator}"


def _monomial(coefficient: Fraction, name: str) -> str:
    """``coefficient * name`` for a non-zero coefficient; ``name`` may be empty."""
    if not name:
        return format_rational(coefficient)
    if coefficient == 1:
        return name
    if coefficient == -1:
        return "-" + name
    return f"{format_rational(coefficient)}*{name}"


def _sum_terms(terms: list[tuple[Fraction, str]]) -> str:
    """``a + b*name - c*name`` for non-zero coefficients, with a common denominator
    pulled out as ``(...)/q`` when every coefficient shares it."""
    if not terms:
        return "0"
    denominators = {c.denominator for c, _ in terms}
    if len(terms) > 1 and denominators != {1} and len(denominators) == 1:
        q = denominators.pop()
        inner = _sum_terms([(Fraction(c.numerator), name) for c, name in terms])
        return f"({inner})/{q}"
    text = _monomial(*terms[0])
    for coefficient, name in terms[1:]:
        if coefficient < 0:
            text += " - " + _monomial(-coefficient, name)
        else:
            text += " + " + _monomial(coefficient, name)
    return text


def _real_basis(F: CyclotomicField) -> list[tuple[str, Value]]:
    """``1`` and the products of ``sqrt2``, ``sqrt3``, ``sqrt5`` that lie in the field."""
    out: list[tuple[str, Value]] = [("", F.one)]
    roots = [
        (name, F.constant(name))
        for name, _, needs in _ROOTS
        if F.conductor % needs == 0
    ]
    for size in (1, 2, 3):
        for chosen in combinations(roots, size):
            name = "*".join(n for n, _ in chosen)
            value = F.one
            for _, v in chosen:
                value = F.mul(value, v)
            out.append((name, value))
    return out


def _imaginary_basis(F: CyclotomicField) -> list[tuple[str, Value]]:
    """``i`` times the real products, when ``i`` itself lies in the field (``4 | M``):
    the grammar embeds ``i`` only then, so a spelling with ``i`` must parse back."""
    if F.conductor % 4:
        return []
    i = F.constant("i")
    return [
        (f"i*{name}" if name else "i", F.mul(i, value))
        for name, value in _real_basis(F)
    ]


def _rational_multiple(
    F: CyclotomicField, x: Value, unit_inverse: Value
) -> Fraction | None:
    """``r`` with ``x = r * unit`` given ``1 / unit``, or ``None``."""
    quotient = F.mul(x, unit_inverse)
    return F.as_rational(quotient) if F.is_rational(quotient) else None


@cache
def _units(conductor: int) -> list[tuple[str, Value]]:
    """The roots that lie in the field with their inverses: ``1 / sqrtN = sqrtN / N``."""
    F = CyclotomicField(conductor)
    out: list[tuple[str, Value]] = [("", F.one)]
    for name, prime, needs in _ROOTS:
        if conductor % needs == 0:
            out.append((name, F.scale(F.constant(name), Fraction(1, prime))))
    return out


@cache
def _nice_basis(conductor: int) -> list[tuple[str, Value]]:
    F = CyclotomicField(conductor)
    return _real_basis(F) + _imaginary_basis(F)


def _coordinates(
    x: Value, basis: list[tuple[str, Value]]
) -> list[tuple[Fraction, str]] | None:
    solution = _solve_rational([value for _, value in basis], x)
    if solution is None:
        return None
    return [(c, name) for c, (name, _) in zip(solution, basis, strict=True) if c]


def pretty(F: CyclotomicField, x: Value) -> str:
    """The shortest spelling among a rational, ``r*zeta^k``, ``r*sqrtN*zeta^k`` and
    ``re + i*im`` over the nice basis; the power-basis polynomial otherwise."""
    if F.is_rational(x):
        return format_rational(x[0])
    candidates: list[str] = []
    M = F.conductor
    for k in range(1, M):
        z_inverse = F.zeta(-k)
        power = "zeta" if k == 1 else f"zeta^{k}"
        for root, unit_inverse in _units(M):
            r = _rational_multiple(F, x, F.mul(unit_inverse, z_inverse))
            if r is not None:
                name = f"{root}*{power}" if root else power
                candidates.append(_monomial(r, name))
    combination = _coordinates(x, _nice_basis(M))
    if combination is not None:
        candidates.append(_sum_terms(combination))
    candidates.append(
        _sum_terms(
            [
                (a, "zeta" if j == 1 else f"zeta^{j}" if j else "")
                for j, a in enumerate(x)
                if a
            ]
        )
    )
    return min(candidates, key=lambda text: (len(text), candidates.index(text)))


# ---- the value grammar: parser -------------------------------------------------------

_TOKEN = re.compile(r"\s*(?:(\d+)|([A-Za-z]\w*)|([-+*/^()]))")
_ATOMS = ("i", "zeta", "sqrt2", "sqrt3", "sqrt5", "phi")


class ValueSyntaxError(ValueError):
    """A string outside the value grammar, or a symbol the conductor cannot hold."""


class _Parser:
    def __init__(self, F: CyclotomicField, text: str) -> None:
        self.F = F
        self.text = text
        self.tokens = self._tokenize(text)
        self.position = 0

    @staticmethod
    def _tokenize(text: str) -> list[tuple[str, str]]:
        tokens: list[tuple[str, str]] = []
        position = 0
        while position < len(text):
            match = _TOKEN.match(text, position)
            if match is None or match.end() == position:
                if text[position:].strip() == "":
                    break
                raise ValueSyntaxError(f"unexpected character in value '{text}'")
            position = match.end()
            integer, name, operator = match.groups()
            if integer is not None:
                tokens.append(("int", integer))
            elif name is not None:
                tokens.append(("name", name))
            elif operator is not None:
                tokens.append(("op", operator))
        return tokens

    def _peek(self) -> tuple[str, str] | None:
        return self.tokens[self.position] if self.position < len(self.tokens) else None

    def _take(self) -> tuple[str, str]:
        token = self._peek()
        if token is None:
            raise ValueSyntaxError(f"unexpected end of value '{self.text}'")
        self.position += 1
        return token

    def parse(self) -> Value:
        value = self._expr()
        if self._peek() is not None:
            raise ValueSyntaxError(f"trailing tokens in value '{self.text}'")
        return value

    def _expr(self) -> Value:
        value = self._term()
        while self._peek() in (("op", "+"), ("op", "-")):
            _, operator = self._take()
            other = self._term()
            value = (
                self.F.add(value, other)
                if operator == "+"
                else self.F.sub(value, other)
            )
        return value

    def _term(self) -> Value:
        value = self._factor()
        while self._peek() in (("op", "*"), ("op", "/")):
            _, operator = self._take()
            other = self._factor()
            if operator == "*":
                value = self.F.mul(value, other)
            else:
                if self.F.is_zero(other):
                    raise ValueSyntaxError(f"division by zero in value '{self.text}'")
                value = self.F.divide(value, other)
        return value

    def _factor(self) -> Value:
        negative = False
        while self._peek() == ("op", "-"):
            self._take()
            negative = not negative
        value = self._atom()
        if self._peek() == ("op", "^"):
            self._take()
            sign = 1
            if self._peek() == ("op", "-"):
                self._take()
                sign = -1
            kind, text = self._take()
            if kind != "int":
                raise ValueSyntaxError(
                    f"an exponent must be an integer in '{self.text}'"
                )
            value = self.F.power(value, sign * int(text))
        return self.F.neg(value) if negative else value

    def _atom(self) -> Value:
        kind, text = self._take()
        if kind == "int":
            return self.F.rational(int(text))
        if kind == "name":
            if text == "zeta":
                return self.F.zeta(1)
            if text in _ATOMS:
                try:
                    return self.F.constant(text)
                except ValueError as exc:
                    raise ValueSyntaxError(str(exc)) from None
            raise ValueSyntaxError(f"unknown symbol '{text}' in value '{self.text}'")
        if (kind, text) == ("op", "("):
            value = self._expr()
            if self._take() != ("op", ")"):
                raise ValueSyntaxError(f"unbalanced parentheses in value '{self.text}'")
            return value
        raise ValueSyntaxError(f"unexpected '{text}' in value '{self.text}'")


def parse(F: CyclotomicField, text: str) -> Value:
    """Parse a string of the value grammar into the field; ``zeta`` is ``zeta_M``."""
    return _Parser(F, text).parse()
