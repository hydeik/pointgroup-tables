"""The ``matrices`` and ``all`` documents (docs/spec/basis-conventions.md
"Verification", docs/spec/output-schema.md, docs/spec/cli.md): ``MX-*``,
``VL-3``, the layout of the matrix table, the union rule of ``all``
(``OS-4``, ``OS-5``), the schema on the seed set, the prototype samples and
the command line."""

from __future__ import annotations

import itertools
import json
import tomllib
from typing import Any

import jsonschema
import pytest
import yaml
from conftest import LARGE_GROUPS, SEED_GROUPS, default_and_large, model
from test_elements import SAMPLES, SCHEMA

from pointgroup_tables import checks
from pointgroup_tables.cli import main
from pointgroup_tables.document import (
    TABLE_ORDER,
    VALUE_FORMS,
    Request,
    build_document,
)
from pointgroup_tables.field import cyclotomic_field
from pointgroup_tables.formats import EMITTERS
from pointgroup_tables.seed import SEED

# the matrix checks that run plainly on every group of default and large
CHECKS = ["MX-1", "MX-3", "MX-4", "MX-5", "MX-7", "MX-8", "VL-3"]
# MX-2 is exhaustive over |G|^2 pairs: slow beyond |G-dagger| = 96
SLOW_PAIRS_ORDER = 48
YAML_LOADER = getattr(yaml, "CSafeLoader", yaml.SafeLoader)


def _mx2_groups() -> list[Any]:
    plain = default_and_large()
    out = []
    for item in plain:
        symbol = item if isinstance(item, str) else item.values[0]
        if symbol in SEED and SEED[symbol].order > SLOW_PAIRS_ORDER:
            out.append(pytest.param(symbol, marks=pytest.mark.slow))
        else:
            out.append(item)
    return out


def _load(fmt: str, text: str) -> Any:
    if fmt == "json":
        return json.loads(text)
    if fmt == "yaml":
        return yaml.load(text, Loader=YAML_LOADER)  # the safe loader, C when built
    return tomllib.loads(text)


# ---- the checks --------------------------------------------------------------------------


@pytest.mark.parametrize("check_id", CHECKS)
@pytest.mark.parametrize("symbol", default_and_large())
def test_matrix_checks(check_id: str, symbol: str) -> None:
    assert checks.check(model(symbol), [check_id]) == []


@pytest.mark.parametrize("symbol", _mx2_groups())
def test_mx2_homomorphism_over_all_pairs(symbol: str) -> None:
    assert checks.check(model(symbol), ["MX-2"]) == []


def test_matrix_checks_apply_to_the_tables_they_need() -> None:
    full = set(checks.applicable_checks(model("C3v")))
    assert {"MX-1", "MX-2", "MX-3", "MX-4", "MX-5", "MX-7", "MX-8", "VL-3"} <= full
    alone = build_document(Request("matrices", "C3v"))
    assert set(checks.applicable_checks(alone)) == {
        "AG-3",
        "CL-2",
        "CL-4",
        "CL-5",
        "DG-6",
        "MX-3",
        "MX-5",
        "MX-7",
        "MX-8",
        "VL-3",
    }
    assert checks.check(alone) == []


# ---- the matrices table -----------------------------------------------------------------


@pytest.mark.parametrize("symbol", SEED_GROUPS)
def test_matrices_document_layout(symbol: str) -> None:
    document = build_document(Request("matrices", symbol))
    assert document["tables"] == ["elements", "matrices"]
    assert document["irreps"] == "all"
    assert list(document)[-2:] == ["elements", "matrices"]
    N = document["group"]["order"]
    characters = build_document(Request("characters", symbol))["characters"]
    entries = document["matrices"]["irreps"]
    assert [entry["index"] for entry in entries] == list(range(len(entries)))
    for entry, record in zip(entries, characters["irreps"], strict=True):
        keys = list(entry)
        assert keys[-2:] == ["exact", "numeric"]
        assert {
            key: value
            for key, value in entry.items()
            if key not in ("exact", "numeric")
        } == record
        d = entry["dimension"]
        for form in ("exact", "numeric"):
            assert len(entry[form]) == N  # over G, also for a projective irrep
            assert all(
                len(matrix) == d and all(len(row) == d for row in matrix)
                for matrix in entry[form]
            )


def test_irreps_filter_keeps_global_indices_in_matrices() -> None:
    projective = build_document(Request("matrices", "C3v", irreps="projective"))
    assert [e["index"] for e in projective["matrices"]["irreps"]] == [3, 4, 5]
    assert projective["irreps"] == "projective"
    ordinary = build_document(Request("matrices", "C3v", irreps="ordinary"))
    assert [e["index"] for e in ordinary["matrices"]["irreps"]] == [0, 1, 2]


@pytest.mark.parametrize(
    "forms", [("exact",), ("numeric",), ("pretty",), ("exact", "numeric", "pretty")]
)
def test_values_option_selects_the_matrix_forms(forms: tuple[str, ...]) -> None:
    document = build_document(Request("matrices", "D3h", values=forms))
    assert document["values"] == list(forms)
    for entry in document["matrices"]["irreps"]:
        present = [form for form in VALUE_FORMS if form in entry]
        assert present == list(forms)
        shapes = {form: [(len(m), len(m[0])) for m in entry[form]] for form in present}
        assert len(set(map(tuple, shapes.values()))) == 1


def test_canonical_generator_images_of_c3v() -> None:
    """The closed forms of the basis-conventions chapter: ``E`` has ``diag(w, w^-1)``
    on ``C3^1`` and ``[[0, 1], [1, 0]]`` on ``sigma_v(1)``; ``E1/2`` has
    ``diag(u, u^-1)`` and ``[[0, i], [i, 0]]`` on the lifted pair (``w = zeta_3``,
    ``u = zeta_6``, ``i = zeta_12^3`` in the document field ``Q(zeta_12)``)."""
    document = build_document(Request("matrices", "C3v", values=("exact",)))
    assert document["field"]["conductor"] == 12
    F = cyclotomic_field(12)
    labels = [item["label"] for item in document["elements"]["items"]]
    a, b = labels.index("C3^1"), labels.index("sigma_v(1)")
    entries = {
        entry["symbol"]: entry["exact"] for entry in document["matrices"]["irreps"]
    }
    zero, one, i = F.exact(F.zero), F.exact(F.one), F.exact(F.zeta(3))
    w, w_inverse = F.exact(F.zeta(4)), F.exact(F.zeta(-4))
    u, u_inverse = F.exact(F.zeta(2)), F.exact(F.zeta(-2))
    assert entries["E"][a] == [[w, zero], [zero, w_inverse]]
    assert entries["E"][b] == [[zero, one], [one, zero]]
    assert entries["E1/2"][a] == [[u, zero], [zero, u_inverse]]
    assert entries["E1/2"][b] == [[zero, i], [i, zero]]
    assert entries["A2"][b] == [[F.exact(F.neg(F.one))]]


# ---- the all document --------------------------------------------------------------------


def _same_options(name: str, request: Request) -> Request:
    """The single-table request under the options ``all`` applied to that table."""
    return Request(
        name,
        request.group,
        values=request.values,
        label_convention=request.label_convention,
        basis_convention=request.basis_convention,
        double=request.double if name in ("elements", "cayley") else False,
        irreps=request.irreps if name in ("characters", "matrices") else "all",
    )


def _assert_union(request: Request, fmt: str) -> None:
    emit = EMITTERS[fmt]
    union = _load(fmt, emit(build_document(request)))
    assert union["tables"] == list(TABLE_ORDER)
    for name in TABLE_ORDER:
        options = _same_options(name, request)
        single = _load(fmt, emit(build_document(options)))
        assert single["tables"] == (
            ["elements"] if name == "elements" else ["elements", name]
        )
        assert single[name] == union[name], (name, fmt)
        if options.double == request.double:
            # a table that takes no --double-group embeds the elements of G
            assert single["elements"] == union["elements"], (name, fmt)
        for key in single:
            if key in ("tables", "irreps", "double") or key in TABLE_ORDER:
                continue
            assert single[key] == union[key], (name, fmt, key)


@pytest.mark.parametrize("fmt", list(EMITTERS))
@pytest.mark.parametrize("symbol", SEED_GROUPS)
def test_os5_all_is_the_union_of_the_five_documents(symbol: str, fmt: str) -> None:
    _assert_union(Request("all", symbol), fmt)


@pytest.mark.parametrize("double", [False, True])
@pytest.mark.parametrize("which", ["ordinary", "projective", "all"])
@pytest.mark.parametrize("symbol", ["C3v", "D3h", "Oh"])
def test_os5_union_holds_under_every_option(
    symbol: str, which: str, double: bool
) -> None:
    _assert_union(
        Request("all", symbol, values=VALUE_FORMS, double=double, irreps=which), "json"
    )


@pytest.mark.parametrize("symbol", SEED_GROUPS)
def test_os4_irrep_records_agree_between_characters_and_matrices(symbol: str) -> None:
    document = build_document(Request("all", symbol))
    records = [
        {key: value for key, value in entry.items() if key not in VALUE_FORMS}
        for entry in document["matrices"]["irreps"]
    ]
    assert records == document["characters"]["irreps"]
    assert list(document)[-5:] == list(TABLE_ORDER)
    assert "irreps" in document and document["irreps"] == "all"


def test_all_double_group_acts_on_elements_and_cayley_only() -> None:
    document = build_document(Request("all", "D3h", double=True))
    N = document["group"]["order"]
    assert document["double"] is True
    assert document["elements"]["group"] == "G_double"
    assert len(document["elements"]["items"]) == 2 * N
    assert document["cayley"]["group"] == "G_double"
    assert len(document["cayley"]["table"]) == 2 * N
    assert document["factor"]["group"] == "G"
    assert len(document["factor"]["table"]) == N
    assert document["characters"]["class_group"] == "G_double"
    assert all(len(entry["exact"]) == N for entry in document["matrices"]["irreps"])
    assert checks.check(document) == []


# ---- the schema (OS-1) on the seed set ---------------------------------------------------


COMBINATIONS = list(
    itertools.product(
        ("ordinary", "projective", "all"),
        (("exact", "numeric"), ("pretty",), ("exact", "numeric", "pretty")),
    )
)


def _os1_cases() -> list[Any]:
    """Every group of seed with every option combination; the icosahedral documents
    (4 MB each) validate plainly under the default options and slowly otherwise."""
    cases = []
    for symbol in SEED_GROUPS:
        for which, forms in COMBINATIONS:
            default = which == "all" and forms == ("exact", "numeric")
            marks = (
                [pytest.mark.slow]
                if SEED[symbol].order > SLOW_PAIRS_ORDER and not default
                else []
            )
            cases.append(pytest.param(symbol, which, forms, marks=marks))
    return cases


@pytest.mark.parametrize(("symbol", "which", "forms"), _os1_cases())
def test_os1_matrices_documents_validate(
    symbol: str, which: str, forms: tuple[str, ...]
) -> None:
    jsonschema.Draft202012Validator(SCHEMA).validate(
        build_document(Request("matrices", symbol, values=forms, irreps=which))
    )


@pytest.mark.parametrize("symbol", SEED_GROUPS)
def test_os1_all_documents_validate(symbol: str) -> None:
    validator = jsonschema.Draft202012Validator(SCHEMA)
    for double in (False, True):
        validator.validate(build_document(Request("all", symbol, double=double)))
    if SEED[symbol].order <= SLOW_PAIRS_ORDER:
        validator.validate(
            build_document(
                Request(
                    "all", symbol, values=("pretty",), double=True, irreps="projective"
                )
            )
        )


@pytest.mark.slow
@pytest.mark.parametrize("symbol", LARGE_GROUPS)
def test_large_matrices_documents_validate(symbol: str) -> None:
    jsonschema.Draft202012Validator(SCHEMA).validate(
        build_document(Request("matrices", symbol))
    )


# ---- the prototype samples ----------------------------------------------------------------


def _assert_like_the_sample(ours: Any, sample: Any, path: str = "") -> None:
    """Equal up to the prototype's numeric noise and its missing ``bethe`` key."""
    if isinstance(sample, dict):
        assert isinstance(ours, dict), path
        assert {k for k in ours if k != "bethe"} == set(sample), path
        for key, value in sample.items():
            _assert_like_the_sample(ours[key], value, f"{path}/{key}")
    elif isinstance(sample, list):
        assert isinstance(ours, list) and len(ours) == len(sample), path
        for k, (a, b) in enumerate(zip(ours, sample, strict=True)):
            _assert_like_the_sample(a, b, f"{path}[{k}]")
    elif isinstance(sample, float):
        assert ours == pytest.approx(sample, abs=1e-15), path
    else:
        assert ours == sample, path


@pytest.mark.parametrize(
    ("sample", "invocation"),
    [
        ("C3v_matrices.json", Request("matrices", "C3v")),
        ("D6h_matrices_exact.json", Request("matrices", "D6h", values=("exact",))),
        ("C3v_all.json", Request("all", "C3v")),
        ("D6h_all.json", Request("all", "D6h")),
    ],
)
def test_matrices_equal_the_prototype_samples(sample: str, invocation: Request) -> None:
    """The matrices table of the samples (the other tables have their own sample
    tests; the Oh and Ih samples stop before the matrices table)."""
    path = SAMPLES / sample
    if not path.is_file():
        pytest.skip(f"the prototype sample {path.name} is not checked in")
    expected = json.loads(path.read_text())
    document = build_document(invocation)
    assert document["tables"] == expected["tables"]
    _assert_like_the_sample(document["matrices"], expected["matrices"], "matrices")
    for key in ("group", "double", "irreps", "values", "field", "basis_convention"):
        assert document[key] == expected[key], key


# ---- the command line ---------------------------------------------------------------------


def test_matrices_and_all_subcommands(
    capsysbinary: pytest.CaptureFixture[bytes],
) -> None:
    assert main(["matrices", "C3v"]) == 0
    document = json.loads(capsysbinary.readouterr().out)
    assert document["tables"] == ["elements", "matrices"]
    assert (
        main(["all", "D3h", "--double-group", "--irreps", "projective", "-f", "toml"])
        == 0
    )
    document = tomllib.loads(capsysbinary.readouterr().out.decode("ascii"))
    assert document["tables"] == list(TABLE_ORDER)
    assert document["double"] is True and document["irreps"] == "projective"
    assert [e["index"] for e in document["matrices"]["irreps"]] == [6, 7, 8]


@pytest.mark.parametrize(
    ("argv", "code"),
    [
        (["matrices", "C3v", "--double-group"], 2),
        (["all", "C3v", "--values", "none"], 2),
        (["matrices", "Xh"], 3),
        (["all", "S3"], 5),
        (["matrices", "C3v", "--basis-convention", "nope"], 4),
    ],
)
def test_matrices_and_all_usage_errors(
    argv: list[str], code: int, capsysbinary: pytest.CaptureFixture[bytes]
) -> None:
    assert main(argv) == code
    assert capsysbinary.readouterr().out == b""


def test_matrices_ih_completes_and_validates() -> None:
    document = build_document(Request("matrices", "Ih"))
    jsonschema.Draft202012Validator(SCHEMA).validate(document)
    assert len(document["matrices"]["irreps"]) == 18
    assert all(len(entry["exact"]) == 120 for entry in document["matrices"]["irreps"])
    assert len(EMITTERS["json"](document).encode("ascii")) > 4_000_000
