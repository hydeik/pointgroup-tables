"""The values layer (docs/spec/values.md): the field, the exact form, the numeric
form against an independent high-precision oracle, and the value grammar."""

from __future__ import annotations

import math
from fractions import Fraction

import mpmath
import pytest
from conftest import SEED_GROUPS, model

from pointgroup_tables import values
from pointgroup_tables.field import (
    CyclotomicField,
    cyclotomic_field,
    cyclotomic_polynomial,
    euler_phi,
)
from pointgroup_tables.values import ValueSyntaxError, numeric, parse, pretty

# ---- the field ------------------------------------------------------------------


@pytest.mark.parametrize(
    ("m", "coefficients"),
    [(1, (-1, 1)), (2, (1, 1)), (3, (1, 1, 1)), (4, (1, 0, 1)), (6, (1, -1, 1)), (12, (1, 0, -1, 0, 1))],
)  # fmt: skip
def test_cyclotomic_polynomials(m: int, coefficients: tuple[int, ...]) -> None:
    assert cyclotomic_polynomial(m) == coefficients
    assert cyclotomic_field(m).degree == euler_phi(m) == len(coefficients) - 1


def test_spec_examples_of_the_exact_form() -> None:
    F = cyclotomic_field(3)
    assert F.exact(F.zeta(1)) == [[0, 1], [1, 1]]
    assert F.exact(F.zeta(-1)) == F.exact(F.zeta(2)) == [[-1, 1], [-1, 1]]
    assert F.exact(F.rational(Fraction(1, 2))) == [[1, 2], [0, 1]]
    assert F.from_exact([[-1, 1], [-1, 1]]) == F.zeta(2)


def test_spec_example_of_the_numeric_form() -> None:
    F = cyclotomic_field(3)
    assert numeric(F, F.zeta(1)) == [-0.5, 0.8660254037844386]
    assert numeric(F, F.zeta(-1)) == [-0.5, -0.8660254037844386]
    assert numeric(F, F.zero) == [0.0, 0.0]


@pytest.mark.parametrize("m", [1, 2, 3, 4, 5, 8, 9, 12, 20, 24, 30])
def test_field_arithmetic(m: int) -> None:
    F = cyclotomic_field(m)
    z = F.zeta(1)
    assert F.power(z, m) == F.one
    assert F.mul(z, F.zeta(-1)) == F.one
    assert F.conjugate(z) == F.zeta(-1)
    for k in range(m):
        x = F.add(F.zeta(k), F.rational(Fraction(1, 3)))
        assert F.mul(x, F.inverse(x)) == F.one
        assert F.divide(x, x) == F.one
        assert F.conjugate(F.conjugate(x)) == x
        assert F.is_real(F.add(x, F.conjugate(x)))
        assert F.is_imaginary(F.sub(x, F.conjugate(x)))
    assert F.power(F.sub(F.zero, F.one), 2) == F.one


def test_roots_of_unity_of_even_order_in_an_odd_conductor() -> None:
    F = cyclotomic_field(5)
    u = F.root_of_unity(10)
    assert F.power(u, 10) == F.one
    assert F.power(u, 5) == F.neg(F.one)
    assert F.power(u, 2) == F.zeta(1)
    assert F.root_of_unity(10, 3) == F.power(u, 3)
    assert F.has_root_of_unity(10) and not F.has_root_of_unity(4)
    with pytest.raises(ValueError):
        F.root_of_unity(4)


def test_named_constants() -> None:
    F = cyclotomic_field(120)
    for name, square in (("sqrt2", 2), ("sqrt3", 3), ("sqrt5", 5)):
        root = F.constant(name)
        assert F.mul(root, root) == F.rational(square)
        assert numeric(F, root)[0] == pytest.approx(math.sqrt(square))
    i = F.constant("i")
    assert F.mul(i, i) == F.neg(F.one)
    phi = F.constant("phi")
    assert F.mul(phi, phi) == F.add(phi, F.one)
    assert numeric(F, phi) == [1.618033988749895, 0.0]
    for conductor, name in ((12, "sqrt2"), (8, "sqrt3"), (12, "sqrt5"), (3, "i")):
        with pytest.raises(ValueError):
            cyclotomic_field(conductor).constant(name)


def test_lift_and_descend() -> None:
    small, large = cyclotomic_field(3), cyclotomic_field(12)
    x = small.add(small.zeta(1), small.rational(2))
    lifted = small.lift(x, large)
    assert lifted == large.add(large.zeta(4), large.rational(2))
    assert large.descend(lifted, small) == x
    assert large.descend(large.zeta(1), small) is None
    with pytest.raises(ValueError):
        large.lift(large.zeta(1), small)


# ---- numeric form against mpmath -------------------------------------------------------


def _mpmath_value(F: CyclotomicField, pairs: list[list[int]]) -> tuple[float, float]:
    with mpmath.workdps(60):
        total = mpmath.mpc(0)
        for j, (p, q) in enumerate(pairs):
            if p:
                total += (
                    mpmath.mpf(p) / q * mpmath.expjpi(mpmath.mpf(2 * j) / F.conductor)
                )
        return float(total.real), float(total.imag)


def _axial_seed() -> list[str]:
    return [
        symbol
        for symbol in SEED_GROUPS
        if symbol not in ("T", "Td", "Th", "O", "Oh", "I", "Ih")
    ]


@pytest.mark.parametrize("symbol", _axial_seed())
def test_numeric_form_agrees_with_mpmath(symbol: str) -> None:
    document = model(symbol)
    F = cyclotomic_field(document["field"]["conductor"])
    table = document["characters"]
    for exact_row, numeric_row in zip(table["exact"], table["numeric"], strict=True):
        for pairs, approx in zip(exact_row, numeric_row, strict=True):
            value = F.from_exact(pairs)
            expected = _mpmath_value(F, pairs)
            if F.is_real(value):
                expected = (expected[0], 0.0)
            if F.is_imaginary(value):
                expected = (0.0, expected[1])
            assert tuple(approx) == expected
            assert math.copysign(1.0, approx[0]) == math.copysign(
                1.0, approx[1]
            ) == 1.0 or (approx[0] != 0 or approx[1] != 0)


def test_numeric_form_never_prints_negative_zero() -> None:
    F = cyclotomic_field(12)
    for k in range(12):
        re, im = numeric(F, F.zeta(k))
        assert (re, math.copysign(1.0, re)) != (0.0, -1.0)
        assert (im, math.copysign(1.0, im)) != (0.0, -1.0)


# ---- the value grammar -------------------------------------------------------------------


@pytest.mark.parametrize(
    ("conductor", "text", "expected"),
    [
        (12, "zeta", "zeta"),
        (12, "zeta^7", "-zeta"),
        (12, "-1 - zeta", "-1 - zeta"),
        (12, "1/2", "1/2"),
        (12, "-2", "-2"),
        (12, "sqrt3", "sqrt3"),
        (12, "zeta + zeta^-1", "sqrt3"),
        (12, "-1/2*i*sqrt3", "-1/2*i*sqrt3"),
        (12, "(1 + sqrt3)/2", "(1 + sqrt3)/2"),
        (12, "i", "i"),
        (3, "zeta^2", "zeta^2"),
        (3, "zeta^4", "zeta"),
        (3, "1/2 + zeta", "1/2 + zeta"),
        (20, "(1 + sqrt5)/4", "(1 + sqrt5)/4"),
        (20, "phi", "(1 + sqrt5)/2"),
        (8, "zeta + zeta^7", "sqrt2"),
        (8, "1/2*sqrt2*zeta", "(1 + i)/2"),
        (8, "sqrt2*zeta^3", "-1 + i"),
    ],
)
def test_pretty_form(conductor: int, text: str, expected: str) -> None:
    F = cyclotomic_field(conductor)
    value = parse(F, text)
    assert pretty(F, value) == expected
    assert parse(F, pretty(F, value)) == value


@pytest.mark.parametrize("conductor", [1, 2, 3, 4, 5, 8, 12, 15, 20, 24, 30])
def test_pretty_form_round_trips_on_every_root_and_sum(conductor: int) -> None:
    F = cyclotomic_field(conductor)
    for k in range(conductor):
        for j in range(conductor):
            value = F.add(F.zeta(k), F.scale(F.zeta(j), Fraction(-1, 2)))
            assert parse(F, pretty(F, value)) == value


@pytest.mark.parametrize(
    ("conductor", "text"),
    [
        (12, "sqrt2"),
        (8, "sqrt3"),
        (3, "i"),
        (3, "i*sqrt3"),
        (12, "zeta +"),
        (12, "foo"),
        (12, "(1"),
        (12, "1/0"),
        (12, "zeta^x"),
        (12, "2 3"),
        (12, ""),
        (12, "1.5"),
    ],
)
def test_parser_rejects_strings_outside_the_grammar(conductor: int, text: str) -> None:
    with pytest.raises(ValueSyntaxError):
        parse(cyclotomic_field(conductor), text)


def test_parser_accepts_the_transcriptions_spellings() -> None:
    F = cyclotomic_field(120)
    assert parse(F, "2*phi - 1") == F.constant("sqrt5")
    assert parse(F, "-(1 - sqrt5)/2") == F.sub(F.constant("phi"), F.one)
    assert parse(F, "sqrt2 * sqrt3") == F.mul(F.constant("sqrt2"), F.constant("sqrt3"))
    assert parse(F, "zeta^120") == F.one
    assert parse(F, "--1") == F.one


def test_format_rational() -> None:
    assert values.format_rational(Fraction(3, 4)) == "3/4"
    assert values.format_rational(Fraction(-2)) == "-2"
