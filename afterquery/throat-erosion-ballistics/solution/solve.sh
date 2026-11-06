#!/bin/bash
set -euo pipefail
mkdir -p /app/output
cat > /tmp/reconstruct.py <<'PY'
import json
import os
import sys
import numpy as np
from scipy.optimize import brentq
from scipy.signal import savgol_filter

DATA = os.environ.get("MOTOR_DATA", "/app/data")
OUT = os.environ.get("MOTOR_OUT", "/app/output/ballistics.json")
MODE = sys.argv[1] if len(sys.argv) > 1 else "full"   # full | const | linear (diagnostic shortcuts)

m = json.load(open(f"{DATA}/motor.json"))
d = np.genfromtxt(f"{DATA}/static_fire.csv", delimiter=",", names=True)
t = d["time_s"]
fs = m["sample_rate_hz"]
pa = m["ambient_pressure_MPa"] * 1e6
g = m["propellant"]["exhaust_gamma"]
rho = m["propellant"]["density_kg_m3"]
m_prop = m["propellant"]["loaded_mass_kg"]
m_ign = m["igniter_gas_mass_kg"]
Do = m["grain"]["outer_diameter_mm"] / 1e3
Di = m["grain"]["initial_bore_diameter_mm"] / 1e3
Lg = m["grain"]["length_mm"] / 1e3
Dt0 = m["nozzle"]["pre_fire_throat_diameter_mm"] / 1e3
Dt1 = m["nozzle"]["post_fire_throat_diameter_mm"] / 1e3
Ae = np.pi / 4 * (m["nozzle"]["exit_diameter_mm"] / 1e3) ** 2
lam = (1 + np.cos(np.deg2rad(m["nozzle"]["divergence_half_angle_deg"]))) / 2
web = (Do - Di) / 2

# 1. sensor corrections: gauge -> absolute, load-cell zero from the pre-fire baseline
p = d["chamber_pressure_MPa_gauge"] * 1e6 + pa
pre = t < -0.05
F = d["thrust_N"] - d["thrust_N"][pre].mean()
p[pre] = pa

# 2. nozzle tables vs throat diameter
def pe_ratio(eps):
    f = lambda M: (1 / M) * ((2 / (g + 1)) * (1 + (g - 1) / 2 * M * M)) ** ((g + 1) / (2 * (g - 1))) - eps
    M = brentq(f, 1.0001, 20.0)
    return (1 + (g - 1) / 2 * M * M) ** (-g / (g - 1))

Dg = np.linspace(0.85 * Dt0, 1.35 * Dt0, 1201)
Ag = np.pi / 4 * Dg**2
per = np.array([pe_ratio(Ae / a) for a in Ag])
cj = np.sqrt(2 * g * g / (g - 1) * (2 / (g + 1)) ** ((g + 1) / (g - 1)) * (1 - per ** ((g - 1) / g)))

# steady burn window: pressure above 40 % of peak
pk = savgol_filter(p, 201, 2)
on = np.flatnonzero(pk > 0.4 * pk.max())
i0, i1 = on[0], on[-1]

def invert_Dt(eta, idx):
    """Throat diameter that reproduces measured thrust at each sample (thrust-coefficient inversion)."""
    model = eta * lam * np.outer(p[idx], Ag * cj) + (np.outer(p[idx], per) - pa) * Ae
    out = np.empty(idx.size)
    for k in range(idx.size):
        out[k] = np.interp(F[idx[k]], model[k], Dg)
    return out

steady = np.arange(i0, i1 + 1)
early = steady[(t[steady] > t[i0] + 0.15) & (t[steady] < t[i0] + 0.55)]  # before the throat insert starts eroding

if MODE == "full":
    # 3. jet efficiency: the pre-erosion throat is known, so calibrate eta there
    eta = brentq(lambda e: np.median(invert_Dt(e, early)) - Dt0, 0.80, 1.05)
    Dt_s = savgol_filter(invert_Dt(eta, steady), int(0.25 * fs) | 1, 2)
else:
    eta = float("nan")
    Dt_s = (np.full(steady.size, Dt0) if MODE == "const"
            else Dt0 + (Dt1 - Dt0) * (t[steady] - t[i0]) / (t[i1] - t[i0]))

Dt = np.full(t.size, Dt0)
Dt[steady] = Dt_s
Dt[i1 + 1:] = Dt_s[-1]
At = np.pi / 4 * Dt**2

# 4. c* from the overall mass balance (propellant + igniter gas leave through the throat).
#    Flow is choked above the critical pressure ratio; below it the subsonic correction applies.
Gam = np.sqrt(g) * (2 / (g + 1)) ** ((g + 1) / (2 * (g - 1)))
crit = (2 / (g + 1)) ** (g / (g - 1))
pr = np.clip(pa / np.maximum(p, pa), 0, 1)
psi = np.where(pr <= crit, 1.0,
               np.sqrt(np.clip(2 * g / (g - 1) * (pr ** (2 / g) - pr ** ((g + 1) / g)), 0, None)) / Gam)
q = p * At * psi                      # = mdot * c*
fire = t > -0.05
cstar = np.trapezoid(q[fire], t[fire]) / (m_prop + m_ign)

# 5. web regression from mass generated (nozzle discharge + gas storage), then r(p)
RT = (Gam * cstar) ** 2
def prop_vol(w):
    Dp, L = Di + 2 * w, Lg - 2 * w
    return np.pi / 4 * (Do**2 - Dp**2) * L
wg = np.linspace(0, web, 4001)
mg = rho * (prop_vol(0) - prop_vol(wg))
out_flow = np.concatenate([[0], np.cumsum(0.5 * (q[1:] + q[:-1]) * np.diff(t))]) / cstar
# propellant burnt = discharged mass + gas stored in the growing free volume - igniter gas
w = np.interp(out_flow - m_ign, mg, wg)
for _ in range(3):
    Vfree = np.pi / 4 * Do**2 * Lg - prop_vol(np.clip(w, 0, web))
    stored = p * Vfree / RT
    w = np.interp(out_flow - m_ign + stored - stored[0], mg, wg)

win = steady[(t[steady] > t[i0] + 0.2) & (w[steady] < 0.97 * web)]
dt = 1 / fs
ws = savgol_filter(w, int(0.2 * fs) | 1, 2)
r = np.gradient(ws, dt)
n, lna = np.polyfit(np.log(p[win] / 4e6), np.log(r[win]), 1)

grid = [1.0, 1.5, 2.0, 2.5, 3.0]
res = {
    "c_star_m_per_s": float(cstar),
    "jet_efficiency": float(eta) if MODE == "full" else 1.0,
    "burn_rate_exponent": float(n),
    "burn_rate_at_4MPa_mm_per_s": float(np.exp(lna) * 1e3),
    "total_impulse_N_s": float(np.trapezoid(F[fire], t[fire])),
    "throat_diameter_mm": {f"{x:.1f}": float(np.interp(x, t, Dt) * 1e3) for x in grid},
}
os.makedirs(os.path.dirname(OUT), exist_ok=True)
json.dump(res, open(OUT, "w"), indent=2)
print(json.dumps(res, indent=2))
PY
python3 /tmp/reconstruct.py
