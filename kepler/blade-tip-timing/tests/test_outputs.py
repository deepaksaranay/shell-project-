import json
from pathlib import Path

import numpy as np
import pytest

OUT = Path("/app/output/blade_modal_results.json")
EXPECTED = json.loads((Path(__file__).parent / "expected.json").read_text())
BLADES = [str(i) for i in range(1, 25)]
F_TOL = 0.0025
Z_TOL = 0.25
Z_MEAN_TOL = 0.07


def _load():
    assert OUT.exists(), f"{OUT} not found"
    return json.loads(OUT.read_text())


def _value(out, blade, key):
    assert "blades" in out and blade in out["blades"], f"blade {blade} missing"
    v = float(out["blades"][blade][key])
    assert np.isfinite(v), f"blade {blade} {key} not finite"
    return v


def test_engine_order():
    out = _load()
    assert int(out["engine_order"]) == EXPECTED["engine_order"]


def test_all_blades_present():
    out = _load()
    missing = [b for b in BLADES if b not in out.get("blades", {})]
    assert not missing, f"missing blades: {missing}"


@pytest.mark.parametrize("blade", BLADES)
def test_resonance_frequency(blade):
    out = _load()
    f = _value(out, blade, "resonance_frequency_Hz")
    ft = EXPECTED["blades"][blade]["resonance_frequency_Hz"]
    assert abs(f / ft - 1) <= F_TOL, f"blade {blade}: {f:.3f} Hz"


@pytest.mark.parametrize("blade", BLADES)
def test_damping_ratio(blade):
    out = _load()
    z = _value(out, blade, "damping_ratio")
    zt = EXPECTED["blades"][blade]["damping_ratio"]
    assert abs(z / zt - 1) <= Z_TOL, f"blade {blade}: {z:.5f}"


def test_damping_fleet_bias():
    out = _load()
    rel = [_value(out, b, "damping_ratio") / EXPECTED["blades"][b]["damping_ratio"] - 1
           for b in BLADES]
    assert abs(np.mean(rel)) <= Z_MEAN_TOL, f"mean relative damping error {np.mean(rel):+.3f}"
