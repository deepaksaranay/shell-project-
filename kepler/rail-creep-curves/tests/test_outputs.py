"""Checks /app/output/creep_forces.json against the exact-theory answer key."""
import json
import math
import os
import stat
from pathlib import Path

import pytest

RESULT = Path("/app/output/creep_forces.json")
KEY = json.loads((Path(__file__).parent / "answer_key.json").read_text())
MU_N = KEY["load_case"]["friction_coefficient"] * KEY["load_case"]["normal_load"]
FORCE_TOL = KEY["limits"]["force_vs_mu_N"] * MU_N
AREA_TOL = KEY["limits"]["area_rel"]
PMAX_TOL = KEY["limits"]["peak_pressure_rel"]
CASE_IDS = sorted(KEY["creep_forces_N"])
TOP_KEYS = {"contact_area_mm2", "max_pressure_MPa", "cases"}
CASE_KEYS = {"id", "Fx_N", "Fy_N"}


def _num(v, what):
    assert isinstance(v, (int, float)) and not isinstance(v, bool), f"{what} is not a number"
    v = float(v)
    assert math.isfinite(v), f"{what} is not finite"
    return v


@pytest.fixture(scope="module")
def result():
    # the result must be a regular file written by the solver, not a link
    st = os.lstat(RESULT) if os.path.lexists(RESULT) else None
    assert st is not None, f"{RESULT} was not written"
    assert stat.S_ISREG(st.st_mode), f"{RESULT} must be a regular file (no symlinks or special files)"
    data = json.loads(RESULT.read_text())
    assert isinstance(data, dict), "top level must be an object"
    assert set(data) == TOP_KEYS, f"top-level keys must be exactly {sorted(TOP_KEYS)}"
    cases = data["cases"]
    assert isinstance(cases, list) and len(cases) == len(CASE_IDS), "cases must list each creepage case once"
    for c in cases:
        assert isinstance(c, dict) and set(c) == CASE_KEYS, f"each case must have exactly {sorted(CASE_KEYS)}"
    return data


def test_contact_area(result):
    a = _num(result["contact_area_mm2"], "contact_area_mm2")
    ref = KEY["normal_contact"]["area_mm2"]
    assert abs(a - ref) <= AREA_TOL * ref, f"contact area {a:.1f} mm2"


def test_max_pressure(result):
    pm = _num(result["max_pressure_MPa"], "max_pressure_MPa")
    ref = KEY["normal_contact"]["peak_pressure_MPa"]
    assert abs(pm - ref) <= PMAX_TOL * ref, f"max pressure {pm:.0f} MPa"


@pytest.mark.parametrize("cid", CASE_IDS)
def test_creep_force(result, cid):
    found = [c for c in result["cases"] if c["id"] == cid]
    assert len(found) == 1, f"case {cid} missing or duplicated"
    fx = _num(found[0]["Fx_N"], f"{cid} Fx_N")
    fy = _num(found[0]["Fy_N"], f"{cid} Fy_N")
    rx, ry = KEY["creep_forces_N"][cid]
    err = math.hypot(fx - rx, fy - ry)
    assert err <= FORCE_TOL, f"case {cid}: force error {err:.0f} N (Fx={fx:.0f}, Fy={fy:.0f})"
