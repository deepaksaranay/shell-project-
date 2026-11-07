"""Checks /app/output/flange_contact.json against the answer key."""
import json
import math
import os
import stat
from pathlib import Path

import pytest

RESULT = Path("/app/output/flange_contact.json")
KEY = json.loads((Path(__file__).parent / "answer_key.json").read_text())
LIM = KEY["limits"]
IDS = sorted(KEY["cases"])
CASE_KEYS = {"id", "lateral_shift_mm", "tread", "flange"}
PATCH_KEYS = {"contact_angle_deg", "normal_force_N", "Fx_N", "Fs_N"}


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
    assert isinstance(data, dict) and set(data) == {"cases"}, "top level must be exactly {'cases': [...]}"
    cases = data["cases"]
    assert isinstance(cases, list) and len(cases) == len(IDS), "one entry per load case"
    out = {}
    for c in cases:
        assert isinstance(c, dict) and set(c) == CASE_KEYS, f"each case must have exactly {sorted(CASE_KEYS)}"
        assert c["id"] in KEY["cases"] and c["id"] not in out, f"unknown or repeated case id {c.get('id')}"
        for k in ("tread", "flange"):
            assert c[k] is None or (isinstance(c[k], dict) and set(c[k]) == PATCH_KEYS), \
                f"case {c['id']} {k} must be null or have exactly {sorted(PATCH_KEYS)}"
        out[c["id"]] = c
    return out


@pytest.mark.parametrize("cid", IDS)
def test_load_case(result, cid):
    c, ref = result[cid], KEY["cases"][cid]
    Q = ref["Q"]
    shift = _num(c["lateral_shift_mm"], f"{cid} lateral_shift_mm")
    assert abs(shift - ref["shift"]) <= LIM["lateral_shift_mm"], f"{cid}: lateral shift {shift:.4f} mm"
    for k in ("tread", "flange"):
        rp, p = ref["patches"][k], c[k]
        if rp is None:
            assert p is None, f"{cid}: the {k} should not be in contact"
            continue
        assert p is not None, f"{cid}: the {k} patch is missing"
        ang, N, Fx, Fs = (_num(p[f], f"{cid} {k} {f}") for f in ("contact_angle_deg", "normal_force_N", "Fx_N", "Fs_N"))
        assert abs(ang - rp[0]) <= LIM["contact_angle_deg"], f"{cid} {k}: contact angle {ang:.2f} deg"
        assert abs(N - rp[1]) <= LIM["normal_force_vs_Q"] * Q, f"{cid} {k}: normal force {N:.0f} N"
        err = math.hypot(Fx - rp[2], Fs - rp[3])
        assert err <= LIM["creep_force_vs_Q"] * Q, f"{cid} {k}: creep force error {err:.0f} N (Fx={Fx:.0f}, Fs={Fs:.0f})"
