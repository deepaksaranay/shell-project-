import json
import math
import os
import pathlib

import pytest

OUT = pathlib.Path(os.environ.get("BALLISTICS_OUTPUT", "/app/output/ballistics.json"))
TRUTH = json.loads((pathlib.Path(__file__).parent / "expected.json").read_text())

SCALAR_KEYS = {
    "c_star_m_per_s",
    "jet_efficiency",
    "burn_rate_exponent",
    "burn_rate_at_4MPa_mm_per_s",
    "total_impulse_N_s",
}
TIME_KEYS = {"1.0", "1.5", "2.0", "2.5", "3.0"}


def _reject_constant(name):
    raise ValueError(f"non-finite JSON constant {name}")


def _is_number(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v)


@pytest.fixture(scope="module")
def result():
    assert OUT.exists(), f"missing {OUT}"
    assert not OUT.is_symlink(), "output must be a regular file, not a symlink"
    assert OUT.is_file()
    data = json.loads(OUT.read_text(), parse_constant=_reject_constant)
    assert isinstance(data, dict), "top level must be a JSON object"
    return data


def test_exact_keys(result):
    assert set(result) == SCALAR_KEYS | {"throat_diameter_mm"}, f"keys: {sorted(result)}"
    th = result["throat_diameter_mm"]
    assert isinstance(th, dict), "throat_diameter_mm must be an object"
    assert set(th) == TIME_KEYS, f"throat_diameter_mm keys: {sorted(th)}"


def test_all_finite_numbers(result):
    for k in SCALAR_KEYS:
        assert _is_number(result[k]), f"{k} is not a finite number"
    for k, v in result["throat_diameter_mm"].items():
        assert _is_number(v), f"throat_diameter_mm[{k}] is not a finite number"


def _rel(key, tol, result):
    got, want = result[key], TRUTH[key]
    assert abs(got - want) <= tol * abs(want), f"{key}: {got:.6g} vs {want:.6g}"


def _abs(key, tol, result):
    got, want = result[key], TRUTH[key]
    assert abs(got - want) <= tol, f"{key}: {got:.6g} vs {want:.6g}"


def test_characteristic_velocity(result):
    _rel("c_star_m_per_s", 0.015, result)


def test_jet_efficiency(result):
    _abs("jet_efficiency", 0.01, result)


def test_burn_rate_exponent(result):
    _abs("burn_rate_exponent", 0.02, result)


def test_burn_rate_at_4mpa(result):
    _rel("burn_rate_at_4MPa_mm_per_s", 0.02, result)


def test_total_impulse(result):
    _rel("total_impulse_N_s", 0.01, result)


@pytest.mark.parametrize("tkey", sorted(TIME_KEYS))
def test_throat_diameter_history(result, tkey):
    got = result["throat_diameter_mm"][tkey]
    want = TRUTH["throat_diameter_mm"][tkey]
    assert abs(got - want) <= 0.05, f"throat diameter at t={tkey} s: {got:.4f} vs {want:.4f} mm"
