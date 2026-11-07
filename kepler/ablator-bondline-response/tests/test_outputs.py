"""Checks /app/output/heatshield_response.json against the answer key."""
import json
import math
import os
import stat
from pathlib import Path

import pytest

RESULT = Path("/app/output/heatshield_response.json")
KEY = json.loads((Path(__file__).parent / "answer_key.json").read_text())
TIMES = KEY["output_times_s"]
QUANTITIES = sorted(KEY["series"])
SNAP_KEYS = {"t_s", *QUANTITIES}
TOP_KEYS = {"snapshots", "peak_bondline_temperature_K"}


def _num(v, what):
    assert isinstance(v, (int, float)) and not isinstance(v, bool), f"{what} is not a number"
    v = float(v)
    assert math.isfinite(v), f"{what} is not finite"
    return v


def _within(value, ref, lim, what):
    tol = lim["abs"] + lim["rel"] * abs(ref)
    assert abs(value - ref) <= tol, f"{what}: {value:.4g} vs expected {ref:.4g} (allowed +-{tol:.3g})"


@pytest.fixture(scope="module")
def result():
    # the result must be a regular file written by the solver, not a link
    st = os.lstat(RESULT) if os.path.lexists(RESULT) else None
    assert st is not None, f"{RESULT} was not written"
    assert stat.S_ISREG(st.st_mode), f"{RESULT} must be a regular file (no symlinks or special files)"
    data = json.loads(RESULT.read_text())
    assert isinstance(data, dict) and set(data) == TOP_KEYS, f"top-level keys must be exactly {sorted(TOP_KEYS)}"
    snaps = data["snapshots"]
    assert isinstance(snaps, list) and len(snaps) == len(TIMES), "one snapshot per output time"
    by_t = {}
    for s in snaps:
        assert isinstance(s, dict) and set(s) == SNAP_KEYS, f"each snapshot must have exactly {sorted(SNAP_KEYS)}"
        t = _num(s["t_s"], "t_s")
        match = [T for T in TIMES if abs(t - T) < 1e-6]
        assert len(match) == 1 and match[0] not in by_t, f"unexpected or repeated output time {t}"
        by_t[match[0]] = s
    return data, by_t


@pytest.mark.parametrize("quantity", QUANTITIES)
def test_time_history(result, quantity):
    _, by_t = result
    lim = KEY["limits"][quantity]
    for T in TIMES:
        v = _num(by_t[T][quantity], f"{quantity} at {T} s")
        _within(v, KEY["series"][quantity][str(T)], lim, f"{quantity} at t = {T} s")


def test_peak_bondline_temperature(result):
    data, _ = result
    v = _num(data["peak_bondline_temperature_K"], "peak_bondline_temperature_K")
    _within(v, KEY["peak_bondline_K"], KEY["limits"]["peak_bondline_K"], "peak bondline temperature")
