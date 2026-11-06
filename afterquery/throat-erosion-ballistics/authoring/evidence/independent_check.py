"""Independent cross-check: forward-model least-squares fit (no thrust-coefficient inversion).

Quasi-steady 0-D ballistics p = (rho*a*Ab*c*/At)^(1/(1-n)) integrated over web, nozzle thrust
model with unknown jet efficiency, and a piecewise-linear throat-diameter history (free knots).
All unknowns (c*, r_4MPa, n, eta, knots) are fitted jointly to measured pressure and thrust,
with the pre/post throat diameters and loaded propellant mass as soft constraints.
Usage: python3 independent_check.py <data_dir> <out.json>
"""
import json
import sys
import numpy as np
from scipy.optimize import brentq, least_squares

D = sys.argv[1]
m = json.load(open(f"{D}/motor.json"))
d = np.genfromtxt(f"{D}/static_fire.csv", delimiter=",", names=True)
t = d["time_s"]
pa = m["ambient_pressure_MPa"] * 1e6
p = d["chamber_pressure_MPa_gauge"] * 1e6 + pa
F = d["thrust_N"] - d["thrust_N"][t < -0.05].mean()
g = m["propellant"]["exhaust_gamma"]
rho = m["propellant"]["density_kg_m3"]
Do, Di, Lg = (m["grain"][k] / 1e3 for k in ("outer_diameter_mm", "initial_bore_diameter_mm", "length_mm"))
Dt0 = m["nozzle"]["pre_fire_throat_diameter_mm"] / 1e3
Dt1 = m["nozzle"]["post_fire_throat_diameter_mm"] / 1e3
Ae = np.pi / 4 * (m["nozzle"]["exit_diameter_mm"] / 1e3) ** 2
lam = (1 + np.cos(np.deg2rad(m["nozzle"]["divergence_half_angle_deg"]))) / 2
web = (Do - Di) / 2

# decimate to 100 Hz, keep the region where the motor is clearly burning
k = slice(None, None, 20)
t, p, F = t[k], p[k], F[k]
burn = p > 0.5 * p.max()
ts, te = t[burn][0] + 0.15, t[burn][-1] - 0.04  # end before burnout; skip the ignition transient (gas storage not modelled)
sel = (t >= ts) & (t <= te)
tt, pp, FF = t[sel], p[sel], F[sel]
knots = np.linspace(ts, te, 15)


def pe_ratio(eps):
    f = lambda M: (1 / M) * ((2 / (g + 1)) * (1 + (g - 1) / 2 * M * M)) ** ((g + 1) / (2 * (g - 1))) - eps
    M = brentq(f, 1.0001, 20.0)
    return (1 + (g - 1) / 2 * M * M) ** (-g / (g - 1))


Dg = np.linspace(0.9 * Dt0, 1.3 * Dt0, 400)
per_g = np.array([pe_ratio(Ae / (np.pi / 4 * x * x)) for x in Dg])


def Ab(w):
    Dp, L = Di + 2 * w, Lg - 2 * w
    return np.pi * Dp * L + np.pi / 2 * (Do**2 - Dp**2)


def model(x):
    cs, r4, n, eta = x[:4]
    Dk = x[4:]
    Dt = np.interp(tt, knots, Dk)
    At = np.pi / 4 * Dt**2
    w = 0.0
    P = np.empty_like(tt)
    # web burnt before the window (ignition rise) from the measured trace
    pre = (t > -0.05) & (t < ts)
    w = np.trapezoid(r4 * 1e-3 * (np.maximum(p[pre], pa) / 4e6) ** n, t[pre])
    w0 = w
    dt = np.diff(tt, prepend=tt[0])
    for i in range(tt.size):
        w += 0 if i == 0 else r4 * 1e-3 * (P[i - 1] / 4e6) ** n * dt[i]
        a = r4 * 1e-3 / (4e6) ** n
        P[i] = (rho * a * Ab(min(w, web)) * cs / At[i]) ** (1 / (1 - n))
    per = np.interp(Dt, Dg, per_g)
    cj = np.sqrt(2 * g * g / (g - 1) * (2 / (g + 1)) ** ((g + 1) / (g - 1)) * (1 - per ** ((g - 1) / g)))
    Fm = eta * lam * P * At * cj + (P * per - pa) * Ae
    # mass balance over the window: discharged gas = propellant consumed between w0 and w_end
    vol = lambda ww: np.pi / 4 * (Do**2 - (Di + 2 * ww) ** 2) * (Lg - 2 * ww)
    m_out = np.trapezoid(P * At / cs, tt)
    m_gen = rho * (vol(w0) - vol(min(w, web)))
    return P, Fm, Dk, (m_out - m_gen, w - web)


def resid(x):
    P, Fm, Dk, (dm, dw) = model(x)
    return np.concatenate([(P - pp) / 2e4, (Fm - FF) / 6.0, [dm / 2e-3, dw / 2e-5],
                           [(Dk[0] - Dt0) / 5e-6, (Dk[-1] - Dt1) / 5e-6],
                           np.diff(Dk, 2) / 2e-5])  # mild smoothness on the knot history


x0 = np.concatenate([[1400, 5.0, 0.3, 0.95], np.linspace(Dt0, Dt1, knots.size)])
sol = least_squares(resid, x0, x_scale=np.concatenate([[10, 0.1, 0.01, 0.01], np.full(knots.size, 1e-5)]))
cs, r4, n, eta = sol.x[:4]
res = {"c_star_m_per_s": cs, "jet_efficiency": eta, "burn_rate_exponent": n,
       "burn_rate_at_4MPa_mm_per_s": r4,
       "throat_diameter_mm": {f"{x:.1f}": float(np.interp(x, knots, sol.x[4:]) * 1e3) for x in [1.0, 1.5, 2.0, 2.5, 3.0]}}
json.dump({k: (float(v) if not isinstance(v, dict) else v) for k, v in res.items()}, open(sys.argv[2], "w"), indent=2)
print(json.dumps({k: (float(v) if not isinstance(v, dict) else v) for k, v in res.items()}, indent=2))
