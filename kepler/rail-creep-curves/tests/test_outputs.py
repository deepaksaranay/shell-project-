"""Checks /app/output/creep_forces.json against the reference rolling-contact solution."""
import json
import math
from pathlib import Path

import pytest

RESULT = Path("/app/output/creep_forces.json")
REF = json.loads((Path(__file__).parent / "reference_solution.json").read_text())
MU_N = REF["mu"] * REF["normal_load_N"]
FORCE_TOL = REF["tolerances"]["force_fraction_of_mu_N"] * MU_N
AREA_TOL = REF["tolerances"]["area_relative"]
PMAX_TOL = REF["tolerances"]["max_pressure_relative"]
CASE_IDS = [c["id"] for c in REF["cases"]]


def _num(v, what):
    assert isinstance(v, (int, float)) and not isinstance(v, bool), f"{what} is not a number"
    v = float(v)
    assert math.isfinite(v), f"{what} is not finite"
    return v


@pytest.fixture(scope="module")
def result():
    assert RESULT.exists(), f"{RESULT} was not written"
    data = json.loads(RESULT.read_text())
    assert isinstance(data, dict), "top level must be an object"
    return data


def test_contact_area(result):
    a = _num(result.get("contact_area_mm2"), "contact_area_mm2")
    ref = REF["contact_area_mm2"]
    assert abs(a - ref) <= AREA_TOL * ref, f"contact area {a:.1f} mm2"


def test_max_pressure(result):
    pm = _num(result.get("max_pressure_MPa"), "max_pressure_MPa")
    ref = REF["max_pressure_MPa"]
    assert abs(pm - ref) <= PMAX_TOL * ref, f"max pressure {pm:.0f} MPa"


@pytest.mark.parametrize("cid", CASE_IDS)
def test_creep_force(result, cid):
    cases = result.get("cases")
    assert isinstance(cases, list), "cases must be a list"
    found = [c for c in cases if isinstance(c, dict) and c.get("id") == cid]
    assert len(found) == 1, f"case {cid} missing or duplicated"
    fx = _num(found[0].get("Fx_N"), f"{cid} Fx_N")
    fy = _num(found[0].get("Fy_N"), f"{cid} Fy_N")
    ref = next(c for c in REF["cases"] if c["id"] == cid)
    err = math.hypot(fx - ref["Fx_N"], fy - ref["Fy_N"])
    assert err <= FORCE_TOL, f"case {cid}: force error {err:.0f} N (Fx={fx:.0f}, Fy={fy:.0f})"
