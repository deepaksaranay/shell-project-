#!/usr/bin/env python3
"""Reference design for the surge-vessel sizing task.

Builds a method-of-characteristics model of the rising main (elastic pipes,
quasi-steady friction, pump rundown from the homologous curves with the rotor
inertia, ideal check valve, polytropic air vessel with liquid level and
differential connection losses), then searches the vessel's initial gas
volume and inflow orifice at a total volume just under the budget, keeping a
safety margin against every limit, and writes the design.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from moc_model import load_system, simulate  # noqa: E402

SYSTEM = "/app/data/system.json"
OUT = "/app/output/surge_vessel_design.json"


def main():
    s = load_system(SYSTEM)
    lim = s["design_limits"]
    V = round(s["design_budget"]["maximum_vessel_total_volume_m3"] - 0.1, 2)
    best = None
    for d in (0.08, 0.12):
        for phi in np.arange(0.38, 0.52, 0.005):
            r = simulate(s, dict(vessel_total_volume_m3=V, initial_gas_volume_m3=V * phi,
                                 inflow_orifice_diameter_m=d), dt=0.01)
            slack = min(r["min_absolute_pressure_head_m"] - lim["minimum_absolute_pressure_head_m"],
                        lim["maximum_gauge_pressure_head_m"] - r["max_gauge_pressure_head_m"],
                        25.0 * (lim["maximum_gas_fraction_of_vessel"] - r["max_gas_fraction"]))
            if best is None or slack > best[0]:
                best = (slack, phi, d, r)
    slack, phi, d, r = best
    design = {"vessel_total_volume_m3": V,
              "initial_gas_volume_m3": round(V * phi, 3),
              "inflow_orifice_diameter_m": d}
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(design, f, indent=2)
    print(json.dumps(design), f"slack={slack:.3f}",
          f"min_abs={r['min_absolute_pressure_head_m']:.2f} max_gauge={r['max_gauge_pressure_head_m']:.2f} "
          f"gas={r['max_gas_fraction']:.3f}")


if __name__ == "__main__":
    main()
