"""The spin double group and the factor system (docs/spec/double-groups.md,
``DG-*`` and the ``G-dagger`` side of ``CL-*``)."""

from __future__ import annotations

from typing import Any

import jsonschema
import pytest
from conftest import (
    DEFAULT_GROUPS,
    LARGE_GROUPS,
    SEED_GROUPS,
    default_and_large,
    double_model,
    model,
)
from test_elements import SCHEMA, load_sample

from pointgroup_tables import checks
from pointgroup_tables.abstract import point_group
from pointgroup_tables.document import Request, build_document
from pointgroup_tables.double import binary_cover, double_group
from pointgroup_tables.seed import SEED
from pointgroup_tables.symbols import parse_group

# ---- the binary rotation groups and the canonical section -------------------


@pytest.mark.parametrize(
    ("symbol", "cover", "order"),
    [("C3", "Z6", 6), ("D4", "Dic4", 16), ("T", "2T", 24), ("O", "2O", 48), ("I", "2I", 120)],
)  # fmt: skip
def test_binary_covers(symbol: str, cover: str, order: int) -> None:
    R = point_group(parse_group(symbol)).rotation_group
    cover_group = binary_cover(R)
    assert cover_group.name == cover
    assert len(cover_group.elements) == order == 2 * len(R.elements)
    ebar = cover_group.ebar
    assert cover_group.mul(ebar, ebar) == cover_group.identity
    for q in cover_group.elements:
        assert cover_group.mul(ebar, q) == cover_group.mul(q, ebar)
    fibres: dict[Any, set[Any]] = {}
    for q in cover_group.elements:
        fibres.setdefault(cover_group.project(q), set()).add(q)
    assert set(fibres) == set(R.elements)
    assert all(len(fibre) == 2 for fibre in fibres.values())
    for r in R.elements:
        lift = cover_group.section(r)
        assert cover_group.project(lift) == r
        assert cover_group.mul(ebar, lift) in fibres[r]
        assert cover_group.mul(ebar, lift) != lift


def test_projection_is_a_homomorphism_onto_the_polyhedral_groups() -> None:
    for symbol in ("T", "O", "I"):
        R = point_group(parse_group(symbol)).rotation_group
        cover = binary_cover(R)
        for p in cover.elements:
            for q in cover.elements:
                assert cover.project(cover.mul(p, q)) == R.mul(
                    cover.project(p), cover.project(q)
                )


def test_canonical_section_lifts_the_positive_generator_and_its_inverse() -> None:
    """``s(C_m^j) = s(C_m)^j`` for ``|j| <= m/2`` and ``w > 0`` on the polyhedral axes."""
    D = double_group(point_group(parse_group("C6")))
    lifted = D.elements[: D.group.order]
    c6 = lifted[1]
    power = D.identity
    for j in range(4):  # E, C6, C3, C2 have angles 0, 60, 120, 180 degrees
        assert lifted[j] == power
        power = D.mul(power, c6)
    assert lifted[4] == D.mul(D.inv(c6), D.inv(c6))  # C3^-1 with angle -120
    assert lifted[5] == D.inv(c6)
    alpha = {"T": 1.0, "O": 2**0.5, "I": (1 + 5**0.5) / 2}
    for symbol in ("T", "O", "I"):
        D = double_group(point_group(parse_group(symbol)))
        for (w, _, _, _), _ in D.elements[1 : D.group.order]:
            real_part = (w[0] + w[1] * alpha[symbol]) / 2  # w = (a + b alpha) / 2
            assert real_part > -1e-12  # every lifted rotation has angle in (-pi, pi]


# ---- the factor table --------------------------------------------------------


def factor_of(symbol: str) -> list[list[int]]:
    return model(symbol)["factor"]["table"]


def test_factor_table_of_c3v_equals_the_prototype_sample() -> None:
    expected = load_sample("C3v_factor.json")
    document = build_document(Request("factor", "C3v"))
    assert document["factor"] == expected["factor"]
    assert document["tables"] == ["elements", "factor"]
    assert document["double"] is False


@pytest.mark.parametrize(
    ("symbol", "sample"),
    [("C3v", "C3v_cayley.json"), ("D6h", "D6h_all.json"), ("Oh", "Oh_all.json"), ("Ih", "Ih_all.json")],
)  # fmt: skip
def test_cayley_and_factor_tables_equal_the_prototype_samples(
    symbol: str, sample: str
) -> None:
    expected = load_sample(sample)
    document = model(symbol)
    assert document["cayley"] == expected["cayley"]
    if "factor" in expected:
        assert document["factor"] == expected["factor"]


def test_binary_rotations_square_to_the_full_turn() -> None:
    items = model("D6h")["elements"]["items"]
    omega = factor_of("D6h")
    for record in items:
        g = record["index"]
        if record["kind"] in ("rotation", "reflection") and record["order"] == 2:
            assert omega[g][g] == -1, record["label"]
    labels = [record["label"] for record in items]
    c3, c3_inverse, c6 = (
        labels.index("C3^1"),
        labels.index("C3^-1"),
        labels.index("C6^1"),
    )
    assert omega[c3][c3] == -1
    assert omega[c3][c3_inverse] == omega[c3_inverse][c3] == 1
    assert omega[c6][c6] == 1
    inversion = labels.index("i")
    assert omega[inversion][inversion] == 1


# ---- DG-* and the G-dagger side of CL-* as pure functions of the model ------


@pytest.mark.parametrize("check_id", ["DG-1", "DG-2", "DG-3", "DG-4", "DG-5"])
@pytest.mark.parametrize("symbol", default_and_large())
def test_factor_system_checks(check_id: str, symbol: str) -> None:
    assert checks.check(model(symbol), [check_id]) == []


@pytest.mark.parametrize(
    "check_id", ["DG-1", "DG-2", "DG-3", "DG-4", "DG-5", "CL-1", "CL-2", "CL-3", "CL-4"]
)  # fmt: skip
@pytest.mark.parametrize("symbol", default_and_large())
def test_double_group_checks(check_id: str, symbol: str) -> None:
    assert checks.check(double_model(symbol), [check_id]) == []


@pytest.mark.parametrize("symbol", SEED_GROUPS)
def test_dg6_seed_table(symbol: str) -> None:
    assert checks.check(model(symbol), ["DG-6"]) == []
    assert checks.check(double_model(symbol), ["DG-6"]) == []
    classes = double_model(symbol)["elements"]["classes"]
    assert len(classes) == SEED[symbol].double_class_count
    assert (
        sum(1 for cls in classes if cls["split"]) == 2 * SEED[symbol].projective_irreps
    )


def test_g_only_checks_skip_a_double_document() -> None:
    document = double_model("C3v")
    assert checks.check(document, ["AG-1", "AG-2", "AG-3", "CL-5"]) == []
    assert checks.check(document) == []


# ---- the double-group element list and classes ---------------------------------


@pytest.mark.parametrize("symbol", default_and_large())
def test_double_group_element_records(symbol: str) -> None:
    base = model(symbol)["elements"]["items"]
    items = double_model(symbol)["elements"]["items"]
    order = len(base)
    assert len(items) == 2 * order
    for record in items:
        i = record["index"]
        original = base[i % order]
        assert record["base_index"] == i % order
        assert record["barred"] is (i >= order)
        assert record["label"] == ("~" if i >= order else "") + original["label"]
        for key in ("parity", "kind", "order", "power", "axis_class", "axis_index"):
            assert record.get(key) == original.get(key)
        assert list(record)[-2:] == ["barred", "base_index"]
    assert items[order]["label"] == "~E"


@pytest.mark.parametrize("symbol", default_and_large())
def test_double_group_classes_are_ordered_by_base_class(symbol: str) -> None:
    base = model(symbol)["elements"]["classes"]
    document = double_model(symbol)
    classes = document["elements"]["classes"]
    order = len(base)
    N = len(model(symbol)["elements"]["items"])
    expected_bases = []
    for cls in base:
        expected_bases.append(cls["index"])
        if cls["regular"]:
            expected_bases.append(cls["index"])
    assert [cls["base_class"] for cls in classes] == expected_bases
    assert len(classes) == order + sum(1 for cls in base if cls["regular"])
    for cls in classes:
        origin = base[cls["base_class"]]
        assert cls["split"] is origin["regular"]
        assert {g % N for g in cls["members"]} == set(origin["members"])
        if cls["split"]:
            assert cls["size"] == origin["size"]
            assert cls["representative"] == origin["representative"] + (
                N if cls["barred"] else 0
            )
            keys = ["index", "label", "size", "representative", "members", "inverse", "base_class", "split", "barred"]  # fmt: skip
        else:
            assert cls["size"] == 2 * origin["size"]
            assert cls["representative"] == origin["representative"]
            keys = ["index", "label", "size", "representative", "members", "inverse", "base_class", "split"]  # fmt: skip
        assert list(cls) == keys
        assert cls["members"][0] == cls["representative"]
        assert cls["members"][1:] == sorted(cls["members"][1:])


@pytest.mark.parametrize(
    ("symbol", "expected"),
    [
        ("C6", "E ~E C6 ~C6 C3 ~C3 C2 ~C2 C3^-1 ~C3^-1 C6^-1 ~C6^-1"),
        ("D3", "E ~E 2C3 2~C3 3C2' 3~C2'"),
        ("D4", "E ~E 2C4 2~C4 C2+~C2 2C2'+2~C2' 2C2''+2~C2''"),
        ("D3h", "E ~E 2C3 2~C3 3C2'+3~C2' 2S3 2~S3 sigma_h+~sigma_h 3sigma_v+3~sigma_v"),
        ("T", "E ~E 3C2+3~C2 4C3 4~C3 4C3^-1 4~C3^-1"),
        ("O", "E ~E 3C2+3~C2 6C2'+6~C2' 8C3 8~C3 6C4 6~C4"),
        ("C3v", "E ~E 2C3 2~C3 3sigma_v 3~sigma_v"),
    ],
)  # fmt: skip
def test_double_group_class_labels(symbol: str, expected: str) -> None:
    classes = double_model(symbol)["elements"]["classes"]
    assert [cls["label"] for cls in classes] == expected.split()


def test_d3_split_classes_mix_lifted_and_barred_elements() -> None:
    """Under the canonical section the class of ``s(C2'(1))`` is ``{C2'(1), ~C2'(2), ~C2'(3)}``."""
    document = double_model("D3")
    labels = [record["label"] for record in document["elements"]["items"]]
    by_label = {cls["label"]: cls for cls in document["elements"]["classes"]}
    assert sorted(labels[g] for g in by_label["3C2'"]["members"]) == [
        "C2'(1)",
        "~C2'(2)",
        "~C2'(3)",
    ]
    assert sorted(labels[g] for g in by_label["3~C2'"]["members"]) == [
        "C2'(2)",
        "C2'(3)",
        "~C2'(1)",
    ]


@pytest.mark.parametrize(
    ("symbol", "sample"),
    [("C3v", "C3v_elements_double.json"), ("D6h", "D6h_elements_double.json"), ("Oh", "Oh_elements_double.json")],
)  # fmt: skip
def test_double_elements_table_equals_the_prototype_sample(
    symbol: str, sample: str
) -> None:
    expected = load_sample(sample)
    document = build_document(Request("elements", symbol, double=True))
    # the prototype read the axis off the label's spelling: every axis of a
    # polyhedral group is a polyhedral axis, and a barred element keeps the axis
    # fields of the element it covers (its `~sigma_h` of D6h came out
    # "polyhedral", its `~sigma_h(1)` of Oh lost the axis index)
    items = expected["elements"]["items"]
    for record in items:
        if "axis_class" in record and symbol == "Oh":
            record["axis_class"] = "polyhedral"
    for record in items:
        if record["barred"]:
            base = items[record["base_index"]]
            for key in ("axis_class", "axis_index"):
                record.pop(key, None)
                if key in base:
                    record[key] = base[key]
            record["inverse"], record["class"] = (
                record.pop("inverse"),
                record.pop("class"),
            )
            for key in ("barred", "base_index"):
                record[key] = record.pop(key)
    assert document["elements"] == expected["elements"]
    assert document["double"] is True
    assert document["tables"] == ["elements"]


# ---- documents ---------------------------------------------------------------


@pytest.mark.parametrize("subcommand", ["cayley", "factor"])
@pytest.mark.parametrize("symbol", SEED_GROUPS)
def test_cayley_and_factor_documents_validate_against_the_schema(
    symbol: str, subcommand: str
) -> None:
    document = build_document(Request(subcommand, symbol))
    jsonschema.Draft202012Validator(SCHEMA).validate(document)
    assert document["tables"] == ["elements", subcommand]
    assert document[subcommand]["group"] == "G"


@pytest.mark.parametrize("subcommand", ["elements", "cayley"])
@pytest.mark.parametrize("symbol", SEED_GROUPS)
def test_double_group_documents_validate_against_the_schema(
    symbol: str, subcommand: str
) -> None:
    document = build_document(Request(subcommand, symbol, double=True))
    jsonschema.Draft202012Validator(SCHEMA).validate(document)
    assert document["double"] is True
    assert document["elements"]["group"] == "G_double"
    if subcommand == "cayley":
        assert document["cayley"]["group"] == "G_double"
        assert len(document["cayley"]["table"]) == 2 * SEED[symbol].order


@pytest.mark.parametrize(
    "symbol",
    [
        *(symbol for symbol in DEFAULT_GROUPS if symbol not in SEED),
        *(pytest.param(symbol, marks=pytest.mark.slow) for symbol in LARGE_GROUPS),
    ],
)
def test_sweep_and_large_double_documents_validate_against_the_schema(
    symbol: str,
) -> None:
    jsonschema.Draft202012Validator(SCHEMA).validate(double_model(symbol))
    jsonschema.Draft202012Validator(SCHEMA).validate(model(symbol))


def test_double_group_is_a_usage_error_on_factor() -> None:
    from pointgroup_tables.errors import UsageError

    with pytest.raises(UsageError):
        build_document(Request("factor", "C3v", double=True))
