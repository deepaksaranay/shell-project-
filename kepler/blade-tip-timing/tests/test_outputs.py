import json
import math
from pathlib import Path

import pytest

OUT = Path("/app/output/blade_modal_results.json")
EXPECTED = json.loads((Path(__file__).parent / "expected.json").read_text())
BLADES = [str(i) for i in range(1, 25)]
F_TOL = 0.0025      # relative, per blade
Z_TOL = 0.25        # relative, per blade
Z_MEAN_TOL = 0.07   # mean signed relative damping error over the 24 blades


def _load():
    assert OUT.exists(), f"{OUT} not found"
    out = json.loads(OUT.read_text())
    assert isinstance(out, dict), "top level must be a JSON object"
    return out


def _value(out, blade, key):
    blades = out.get("blades")
    assert isinstance(blades, dict) and blade in blades, f"blade {blade} missing"
    entry = blades[blade]
    assert isinstance(entry, dict) and key in entry, f"blade {blade}: {key} missing"
    v = entry[key]
    assert isinstance(v, (int, float)) and not isinstance(v, bool), f"blade {blade}: {key} not a number"
    v = float(v)
    assert math.isfinite(v), f"blade {blade}: {key} not finite"
    return v


def test_engine_order():
    out = _load()
    eo = out.get("engine_order")
    assert isinstance(eo, int) and not isinstance(eo, bool), "engine_order must be an integer"
    assert eo == EXPECTED["engine_order"]


def test_all_blades_present():
    out = _load()
    blades = out.get("blades")
    assert isinstance(blades, dict), "blades must be an object"
    missing = [b for b in BLADES if b not in blades]
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
    assert z > 0, f"blade {blade}: damping must be positive"
    assert abs(z / zt - 1) <= Z_TOL, f"blade {blade}: {z:.5f}"


def test_damping_fleet_bias():
    out = _load()
    rel = [_value(out, b, "damping_ratio") / EXPECTED["blades"][b]["damping_ratio"] - 1
           for b in BLADES]
    mean = sum(rel) / len(rel)
    assert abs(mean) <= Z_MEAN_TOL, f"mean relative damping error {mean:+.3f}"
