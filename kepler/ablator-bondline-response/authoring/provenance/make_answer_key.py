#!/usr/bin/env python3
"""Write tests/answer_key.json from the converged reference run (Na = 480)."""
import json, os
here = os.path.dirname(os.path.abspath(__file__))
ref = json.load(open(os.path.join(here, "truth_Na480.json")))
q = ["recession_mm", "surface_temperature_K", "bondline_temperature_K",
     "pyrolysis_front_depth_mm", "char_front_depth_mm", "mass_loss_kg_per_m2"]
key = {
    "output_times_s": [s["t_s"] for s in ref["snapshots"]],
    "series": {k: {str(s["t_s"]): round(s[k], 6) for s in ref["snapshots"]} for k in q},
    "peak_bondline_K": round(ref["peak"]["peak_bondline_temperature_K"], 4),
    "limits": {
        "recession_mm": {"rel": 0.03, "abs": 0.01},
        "surface_temperature_K": {"rel": 0.0, "abs": 12.0},
        "bondline_temperature_K": {"rel": 0.0, "abs": 3.0},
        "pyrolysis_front_depth_mm": {"rel": 0.03, "abs": 0.1},
        "char_front_depth_mm": {"rel": 0.03, "abs": 0.1},
        "mass_loss_kg_per_m2": {"rel": 0.02, "abs": 0.01},
        "peak_bondline_K": {"rel": 0.0, "abs": 3.0}
    },
    "provenance": "authoring/provenance/ablation.py, contracting-grid finite volumes, 480 ablator nodes, BDF rtol 1e-7 (240 nodes differ by < 0.3 K and < 0.3%)"
}
json.dump(key, open(os.path.join(here, "..", "..", "tests", "answer_key.json"), "w"), indent=1)
print("ok")
