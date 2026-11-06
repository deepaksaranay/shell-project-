import json, sys, numpy as np
sys.path.insert(0, "bearing-impedance-identification/solution")
from solve import identify
truth = json.load(open("bearing-impedance-identification/tests/expected.json"))
D = "bearing-impedance-identification/environment/data"
variants = {"reference": {}, "rect_window": {"rect": True}, "keep_all_tones": {"keep_all": True},
  "rect+keep_all": {"rect": True, "keep_all": True},
  "no_baseline": {"no_baseline": True}, "relative_for_support": {"relative_for_support": True},
  "no_rotation": {"no_rotation": True}, "no_probe_sign": {"no_probe_sign": True},
  "diag_force": {"diag_force": True}, "no_inertia": {"no_inertia": True}, "no_mass": {"no_mass": True}}
for name, o in variants.items():
    r = identify(D, o)
    worstK = worstC = 0
    for c, v in truth.items():
        Kt, Ct = np.array(v["K"]), np.array(v["C"])
        K, C = np.array(r[c]["K"]), np.array(r[c]["C"])
        sK = max(abs(Kt[0,0]), abs(Kt[1,1])); sC = max(abs(Ct[0,0]), abs(Ct[1,1]))
        worstK = max(worstK, np.max(np.abs(K-Kt))/sK); worstC = max(worstC, np.max(np.abs(C-Ct))/sC)
    print(f"{name:22s} worst |dK|/Kdiag={worstK:.4f}  worst |dC|/Cdiag={worstC:.4f}")
