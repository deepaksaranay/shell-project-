import json
import math
from pathlib import Path

import pytest

from hydraulics import load_system, simulate

OUT = Path("/app/output/surge_vessel_design.json")
SYSTEM = load_system(Path(__file__).parent / "system.json")
LIM = SYSTEM["design_limits"]
BUDGET = SYSTEM["design_budget"]["maximum_vessel_total_volume_m3"]
KEYS = ("vessel_total_volume_m3", "initial_gas_volume_m3", "inflow_orifice_diameter_m")
DT = 0.005


def _design():
    assert OUT.exists(), f"{OUT} not found"
    d = json.loads(OUT.read_text())
    assert isinstance(d, dict), "design must be a JSON object"
    out = {}
    for k in KEYS:
        v = d.get(k)
        assert isinstance(v, (int, float)) and not isinstance(v, bool), f"{k} must be a number"
        v = float(v)
        assert math.isfinite(v), f"{k} must be finite"
        out[k] = v
    return out


@pytest.fixture(scope="session")
def result():
    d = _design()
    return d, simulate(SYSTEM, d, dt=DT)


def test_design_within_budget_and_ranges():
    d = _design()
    assert 0 < d["vessel_total_volume_m3"] <= BUDGET, f"vessel volume {d['vessel_total_volume_m3']} m3"
    assert 0 < d["initial_gas_volume_m3"] <= LIM["maximum_gas_fraction_of_vessel"] * d["vessel_total_volume_m3"], \
        "initial gas volume out of range"
    lo, hi = LIM["inflow_orifice_diameter_range_m"]
    assert lo <= d["inflow_orifice_diameter_m"] <= hi, "inflow orifice diameter out of range"


def test_minimum_pressure(result):
    d, r = result
    assert r["min_absolute_pressure_head_m"] >= LIM["minimum_absolute_pressure_head_m"], \
        f"minimum absolute pressure head {r['min_absolute_pressure_head_m']:.2f} m at (t, chainage) {r['min_pressure_time_chainage']}"


def test_maximum_pressure(result):
    d, r = result
    assert r["max_gauge_pressure_head_m"] <= LIM["maximum_gauge_pressure_head_m"], \
        f"maximum gauge pressure head {r['max_gauge_pressure_head_m']:.2f} m at (t, chainage) {r['max_pressure_time_chainage']}"


def test_vessel_gas_fraction(result):
    d, r = result
    assert r["max_gas_fraction"] <= LIM["maximum_gas_fraction_of_vessel"], \
        f"gas fills {r['max_gas_fraction']:.3f} of the vessel"
