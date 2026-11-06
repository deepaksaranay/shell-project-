# Run from the bundle root: refine the minimum feasible volume (verifier model, dt=0.005).
import sys, numpy as np, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, "tests"); import hydraulics as h
s = h.load_system("environment/data/system.json")
DT = float(sys.argv[1]) if len(sys.argv) > 1 else 0.005
def ev(V, phi, d):
    return h.simulate(s, dict(vessel_total_volume_m3=V, initial_gas_volume_m3=V*phi, inflow_orifice_diameter_m=d), dt=DT)
for V in [float(v) for v in sys.argv[2].split(",")]:
    best = None
    for phi in np.arange(0.40, 0.51, 0.01):
        r = ev(V, phi, 0.08)
        ok = h.feasible(r, s)
        sl = min(r["min_absolute_pressure_head_m"] - 3.0, 100 * (0.9 - r["max_gas_fraction"]), 105 - r["max_gauge_pressure_head_m"])
        if best is None or sl > best[0]:
            best = (sl, phi, r, ok)
    sl, phi, r, ok = best
    print(f"V={V:.2f} phi*={phi:.2f} slack={sl:+.3f} feasible={ok} min_abs={r['min_absolute_pressure_head_m']:.3f} gas={r['max_gas_fraction']:.4f} maxg={r['max_gauge_pressure_head_m']:.2f}", flush=True)
