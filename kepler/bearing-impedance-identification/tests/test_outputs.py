import json
from pathlib import Path

import numpy as np
import pytest

OUT = Path("/app/output/bearing_coefficients.json")
EXPECTED = json.loads((Path(__file__).parent / "expected.json").read_text())
TOL = 0.03
CONDS = sorted(EXPECTED)


def _load():
    assert OUT.exists(), f"{OUT} not found"
    return json.loads(OUT.read_text())


def _matrix(entry, key):
    assert key in entry, f"missing {key}"
    m = np.asarray(entry[key], dtype=float)
    assert m.shape == (2, 2), f"{key} must be 2x2, got {m.shape}"
    assert np.all(np.isfinite(m)), f"{key} has non-finite values"
    return m


def test_all_conditions_present():
    out = _load()
    missing = [c for c in CONDS if c not in out]
    assert not missing, f"missing conditions: {missing}"


@pytest.mark.parametrize("cond", CONDS)
def test_stiffness(cond):
    out = _load()
    assert cond in out
    K = _matrix(out[cond], "K_N_per_m")
    Kt = np.asarray(EXPECTED[cond]["K"])
    scale = max(abs(Kt[0, 0]), abs(Kt[1, 1]))
    err = np.abs(K - Kt) / scale
    assert err.max() <= TOL, f"{cond} K error/scale = {err.round(4).tolist()}"


@pytest.mark.parametrize("cond", CONDS)
def test_damping(cond):
    out = _load()
    assert cond in out
    C = _matrix(out[cond], "C_Ns_per_m")
    Ct = np.asarray(EXPECTED[cond]["C"])
    scale = max(abs(Ct[0, 0]), abs(Ct[1, 1]))
    err = np.abs(C - Ct) / scale
    assert err.max() <= TOL, f"{cond} C error/scale = {err.round(4).tolist()}"
