# Run from the bundle root. Finds the smallest feasible vessel under the verifier model.
import sys, json, time, numpy as np
sys.path.insert(0, "tests"); import hydraulics as h
s = h.load_system("environment/data/system.json")
DT = float(sys.argv[1]) if len(sys.argv) > 1 else 0.02

def ev(V, phi, d, dt=DT):
    return h.simulate(s, dict(vessel_total_volume_m3=V, initial_gas_volume_m3=V*phi, inflow_orifice_diameter_m=d), dt=dt)

def best_at(V):
    """Return the most feasible (phi, d) at volume V and its slack."""
    best = None
    for phi in np.linspace(0.30, 0.80, 11):
        for d in [0.06, 0.08, 0.10, 0.13, 0.16, 0.20, 0.25, 0.35]:
            r = ev(V, phi, d)
            lim = s["design_limits"]
            slack = min(r["min_absolute_pressure_head_m"] - lim["minimum_absolute_pressure_head_m"],
                        lim["maximum_gauge_pressure_head_m"] - r["max_gauge_pressure_head_m"],
                        100 * (lim["maximum_gas_fraction_of_vessel"] - r["max_gas_fraction"]))
            if best is None or slack > best[0]:
                best = (slack, phi, d, r)
    return best

t = time.time()
for V in [float(v) for v in (sys.argv[2] if len(sys.argv) > 2 else "10,12,14,16,18,20").split(",")]:
    sl, phi, d, r = best_at(V)
    print(f"V={V} best slack={sl:.2f} phi={phi:.2f} d={d} min_abs={r['min_absolute_pressure_head_m']:.2f} max_g={r['max_gauge_pressure_head_m']:.2f} gas={r['max_gas_fraction']:.3f}  ({time.time()-t:.0f}s)", flush=True)
