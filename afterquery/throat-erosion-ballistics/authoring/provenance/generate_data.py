"""Deterministic generator for the throat-erosion-ballistics task.

0-D (lumped) internal ballistics of a single-segment BATES solid motor with
chamber gas storage, pressure-driven throat erosion after a thermal onset
delay, a supersonic conical nozzle with divergence and jet-efficiency losses,
and realistic sensor errors (gauge pressure, load-cell zero offset, noise).
Writes the agent-visible data to environment/data and the truth to
tests/expected.json.
"""
import json
import pathlib
import numpy as np
from scipy.optimize import brentq

ROOT = pathlib.Path(__file__).resolve().parents[2]
DATA = ROOT / "environment" / "data"

# ---- propellant / gas (truth) ----
rho_p = 1765.0           # kg/m^3 (given)
gamma = 1.21             # (given)
c_star = 1492.0          # m/s, delivered (unknown to solver)
n_true = 0.36            # burn-rate exponent (unknown)
r_ref = 7.9e-3           # m/s at p_ref (unknown)
p_ref = 7.0e6
# ---- grain (given) ----
D_o, D_i, L_g = 0.098, 0.040, 0.300
web = (D_o - D_i) / 2
V_case = np.pi / 4 * D_o**2 * (L_g + 0.030)
# ---- nozzle ----
D_t0 = 0.01800           # m (given, pre-fire)
D_e = 0.04600            # m (given, not eroding)
alpha = np.deg2rad(15.0)
lam = (1 + np.cos(alpha)) / 2
eta = 0.952              # jet efficiency (unknown)
# ---- erosion (unknown) ----
t_on = 0.85              # s, throat insert heat-up before erosion starts
e_ref = 0.55e-3          # m/s diameter growth at p_ref
e_exp = 0.85
p_amb = 0.0942e6         # Pa (given)

Gam = np.sqrt(gamma) * (2 / (gamma + 1)) ** ((gamma + 1) / (2 * (gamma - 1)))
RT = (Gam * c_star) ** 2


def nozzle_flow(p, At):
    """Mass flow through the throat: choked (p*At/c*) or subsonic when pa/p exceeds the critical ratio."""
    pr = p_amb / p
    crit = (2 / (gamma + 1)) ** (gamma / (gamma - 1))
    if pr <= crit:
        return p * At / c_star
    psi = np.sqrt(2 * gamma / (gamma - 1) * (pr ** (2 / gamma) - pr ** ((gamma + 1) / gamma)))
    return p * At / c_star * psi / Gam


def burn_area(w):
    if w >= web:
        return 0.0
    Dp, L = D_i + 2 * w, L_g - 2 * w
    return np.pi * Dp * L + 2 * np.pi / 4 * (D_o**2 - Dp**2)


def v_free(w):
    w = min(w, web)
    Dp, L = D_i + 2 * w, L_g - 2 * w
    return V_case - np.pi / 4 * (D_o**2 - Dp**2) * L


def area_ratio_to_pe(eps):
    g = gamma

    def f(M):
        return (1 / M) * ((2 / (g + 1)) * (1 + (g - 1) / 2 * M * M)) ** ((g + 1) / (2 * (g - 1))) - eps
    M = brentq(f, 1.0001, 20.0)
    return (1 + (g - 1) / 2 * M * M) ** (-g / (g - 1))


def thrust(p, Dt):
    if p <= p_amb * 1.05:
        return 0.0
    At, Ae = np.pi / 4 * Dt**2, np.pi / 4 * D_e**2
    pe = p * area_ratio_to_pe(Ae / At)
    jet = p * At * np.sqrt(2 * gamma**2 / (gamma - 1) * (2 / (gamma + 1)) ** ((gamma + 1) / (gamma - 1))
                           * (1 - (pe / p) ** ((gamma - 1) / gamma)))
    return max(0.0, eta * lam * jet + (pe - p_amb) * Ae)


def simulate(dt=2e-5, t_end=6.0):
    p, w, Dt, t = p_amb, 0.0, D_t0, 0.0
    out = []
    m_burnt = 0.0
    while t < t_end:
        Ab = burn_area(w)
        r = r_ref * (p / p_ref) ** n_true if Ab > 0 else 0.0
        mdot_ign = 0.25 if t < 0.04 else 0.0
        At = np.pi / 4 * Dt**2
        mgen = rho_p * Ab * r
        V = v_free(w)
        dp = (mgen + mdot_ign - nozzle_flow(p, At) - p / RT * Ab * r) * RT / V
        dDt = e_ref * (p / p_ref) ** e_exp if (t > t_on and Ab > 0) else 0.0
        out.append((t, p, w, Dt, Ab, r))
        p = max(p_amb, p + dp * dt)
        w += r * dt
        Dt += dDt * dt
        m_burnt += mgen * dt
        t += dt
    return np.array(out), m_burnt


def main():
    rng = np.random.default_rng(1113)
    DATA.mkdir(parents=True, exist_ok=True)
    sim, m_burnt = simulate()
    fs = 2000.0
    t_rec = np.arange(-0.5, 6.0, 1 / fs)
    ts = sim[:, 0]
    p_true = np.interp(t_rec, ts, sim[:, 1], left=p_amb)
    Dt_true = np.interp(t_rec, ts, sim[:, 3], left=D_t0)
    F_true = np.array([thrust(p, d) for p, d in zip(p_true, Dt_true)])
    p_gauge_MPa = (p_true - p_amb) / 1e6 + rng.normal(0, 0.012, t_rec.size)
    F_meas = F_true + 41.7 + rng.normal(0, 6.0, t_rec.size)  # load-cell zero offset + noise
    np.savetxt(DATA / "static_fire.csv", np.column_stack([t_rec, p_gauge_MPa, F_meas]),
               delimiter=",", header="time_s,chamber_pressure_MPa_gauge,thrust_N", comments="",
               fmt=["%.4f", "%.5f", "%.2f"])
    m_prop = rho_p * np.pi / 4 * (D_o**2 - D_i**2) * L_g
    rig = {
        "sample_rate_hz": fs,
        "ambient_pressure_MPa": p_amb / 1e6,
        "propellant": {"density_kg_m3": rho_p, "loaded_mass_kg": round(m_prop, 4),
                        "exhaust_gamma": gamma},
        "grain": {"type": "single cylindrical BATES segment: bore and both end faces burn, outer surface bonded/inhibited",
                  "outer_diameter_mm": D_o * 1e3, "initial_bore_diameter_mm": D_i * 1e3,
                  "length_mm": L_g * 1e3},
        "igniter_gas_mass_kg": 0.010,
        "nozzle": {"pre_fire_throat_diameter_mm": D_t0 * 1e3,
                   "post_fire_throat_diameter_mm": round(sim[-1, 3] * 1e3, 3),
                   "exit_diameter_mm": D_e * 1e3, "divergence_half_angle_deg": 15.0},
        "channels": {"chamber_pressure_MPa_gauge": "chamber pressure, gauge (relative to ambient), MPa",
                     "thrust_N": "load-cell thrust, N"},
    }
    (DATA / "motor.json").write_text(json.dumps(rig, indent=2) + "\n")
    # truth at the graded instants
    grid = [1.0, 1.5, 2.0, 2.5, 3.0]
    truth = {
        "c_star_m_per_s": c_star,
        "jet_efficiency": eta,
        "burn_rate_exponent": n_true,
        "burn_rate_at_4MPa_mm_per_s": r_ref * 1e3 * (4.0e6 / p_ref) ** n_true,
        "total_impulse_N_s": float(np.trapezoid(F_true, t_rec)),
        "throat_diameter_mm": {f"{t:.1f}": float(np.interp(t, ts, sim[:, 3]) * 1e3) for t in grid},
        "burnout_time_s": float(ts[np.argmax(sim[:, 2] >= web)]),
        "mass_burnt_kg": m_burnt,
    }
    (ROOT / "tests" / "expected.json").write_text(json.dumps(truth, indent=2) + "\n")
    print(json.dumps(truth, indent=2))
    print("peak p MPa", sim[:, 1].max() / 1e6, "peak F", F_true.max())


if __name__ == "__main__":
    main()
