#!/usr/bin/env python3
"""Generate the example dataset (shipped) and the held-out datasets + truth
(sealed in tests/) for the blade-tip-timing task.

Lives in authoring/ (never mounted). Run from the bundle root:
    python3 authoring/provenance/generate.py            # all configs
    BTT_ROOT=/tmp/x BTT_SEED_OFFSET=17 python3 ...      # re-seeded copies

Model (per dataset): each blade i is an uncoupled SDOF oscillator
    x'' + 2 z_i wn_i x' + wn_i^2 x = wn_i^2 a_i cos(EO*Phi(t) - EO*2pi*i/B + psi)
    fn_i(Omega)^2 = f0_i^2 + S Omega^2          (Southwell, Omega in rev/s)
integrated (RK4) through the speed history; tip arrival times at casing
probes solve theta_i(t) + x_i(t)/R = alpha_p. Added: probe-angle install
errors, stagger, untwist, TOA jitter, missed/spurious probe triggers, a probe
dropout window, missed/spurious OPR pulses.

Truth = generator parameters (independent of the reference solution):
EO, f_res_i = EO f0_i / sqrt(EO^2 - S), zeta_i, a_i / (2 zeta_i).
"""
import json
import os

import numpy as np

ROOT = os.environ.get("BTT_ROOT", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
SEED_OFFSET = int(os.environ.get("BTT_SEED_OFFSET", "0"))
DT = 1.0 / 25000.0
JIT_TOA = 0.30e-6
JIT_OPR = 0.15e-6

# name: (seed, blades, EO, Southwell S, tip radius m, n probes, rpm start, rpm end,
#        duration s, fixed probe angles or None, dropout probe index (1-based))
CONFIGS = {
    "example": (7710, 24, 6, 4.2, 0.2475, 5, 6550.0, 7550.0, 6.0, [0.0, 14.0, 48.5, 56.5, 153.0], 4),
    "rotor_a": (1101, 22, 7, 3.6, 0.2310, 6, 5400.0, 6400.0, 5.0, None, 2),
    "rotor_b": (2202, 26, 5, 4.8, 0.2620, 5, 7600.0, 8800.0, 7.0, None, 5),
    "rotor_c": (3303, 24, 8, 3.0, 0.2475, 6, 7200.0, 6200.0, 6.0, None, 1),  # run-down
}


def layout_quality(al_deg, eo, eo_max=30):
    a = np.deg2rad(al_deg)
    A = lambda e: np.column_stack([np.ones(len(a)), np.cos(e * a), np.sin(e * a)])
    P = A(eo)
    noise = np.sqrt(np.diag(np.linalg.inv(P.T @ P))[1:].sum())
    worst = 1.0
    for e in range(1, eo_max + 1):
        if e == eo:
            continue
        Pe = A(e)
        proj = Pe @ np.linalg.pinv(Pe)
        for ph in np.linspace(0, np.pi, 12):
            y = np.cos(eo * a + ph)
            r = y - proj @ y
            yc = y - y.mean()
            worst = min(worst, (r @ r) / (yc @ yc))
    return noise, worst


def choose_layout(rng, n, eo):
    best = None
    for _ in range(1500):
        al = np.sort(np.round(rng.uniform(0, 200, n), 1))
        al -= al[0]
        if np.min(np.diff(al)) < 6.0:
            continue
        noise, worst = layout_quality(al, eo)
        sc = noise / np.sqrt(worst)
        if best is None or sc < best[0]:
            best = (sc, al)
    return best[1]


def generate(name, data_dir, truth_path):
    seed, B, EO, S, R, NPR, RPM0, RPM1, TDUR, angles, gap_probe = CONFIGS[name]
    seed += SEED_OFFSET
    rng = np.random.default_rng(seed)
    rng_lay = np.random.default_rng(seed + 7)
    probe_nom = np.array(angles if angles is not None else choose_layout(rng_lay, NPR, EO), float)
    probe_err = rng.uniform(-0.08, 0.08, NPR)
    om_mid = 0.5 * (RPM0 + RPM1) / 60.0
    f0 = om_mid * np.sqrt(EO**2 - S) * (1 + 0.012 * rng.standard_normal(B))
    zeta = rng.uniform(0.0012, 0.0035, B)
    amp = 1.5e-6 * (1 + 0.15 * rng.standard_normal(B))
    stagger = 0.4e-3 * rng.standard_normal(B)
    untwist = 0.25e-3 * (1 + 0.3 * rng.standard_normal(B))
    psi = rng.uniform(0, 2 * np.pi)
    w1, w2 = rng.uniform(0.25, 0.5), rng.uniform(1.5, 2.3)

    def speed_hz(t):
        lin = (RPM0 + (RPM1 - RPM0) * t / TDUR) / 60.0
        return lin + 0.04 * np.sin(2 * np.pi * w1 * t) + 0.02 * np.sin(2 * np.pi * w2 * t + 1.0)

    nt = int(TDUR / DT) + 1
    t = np.arange(nt) * DT
    Om = speed_hz(t)
    Phi = np.concatenate([[0.0], np.cumsum(0.5 * (Om[1:] + Om[:-1]) * DT)]) * 2 * np.pi + 0.3

    blade_ph = EO * 2 * np.pi * np.arange(B) / B
    w0sq = (2 * np.pi * f0) ** 2

    def rhs(Phi_v, Om_v, x, v):
        wn2 = w0sq + S * (2 * np.pi * Om_v) ** 2
        f = wn2 * amp * np.cos(EO * Phi_v - blade_ph + psi)
        return v, f - 2 * zeta * np.sqrt(wn2) * v - wn2 * x

    X = np.zeros((nt, B)); V = np.zeros((nt, B))
    x = np.zeros(B); v = np.zeros(B)
    for k in range(nt - 1):
        p0, p1 = Phi[k], Phi[k + 1]; o0, o1 = Om[k], Om[k + 1]
        pm, om = 0.5 * (p0 + p1), 0.5 * (o0 + o1)
        k1x, k1v = rhs(p0, o0, x, v)
        k2x, k2v = rhs(pm, om, x + 0.5 * DT * k1x, v + 0.5 * DT * k1v)
        k3x, k3v = rhs(pm, om, x + 0.5 * DT * k2x, v + 0.5 * DT * k2v)
        k4x, k4v = rhs(p1, o1, x + DT * k3x, v + DT * k3v)
        x = x + DT / 6 * (k1x + 2 * k2x + 2 * k3x + k4x)
        v = v + DT / 6 * (k1v + 2 * k2v + 2 * k3v + k4v)
        X[k + 1], V[k + 1] = x, v

    os.makedirs(data_dir, exist_ok=True)
    # OPR with faults
    n_rev = int(np.floor((Phi[-1] - 0.05) / (2 * np.pi)))
    opr_t = np.interp(2 * np.pi * np.arange(1, n_rev + 1), Phi, t)
    opr_t += JIT_OPR * rng.standard_normal(len(opr_t))
    rng_opr = np.random.default_rng(seed + 1)
    n = len(opr_t)
    drop = rng_opr.choice(np.arange(10, n - 10), 4, replace=False)
    spur_i = rng_opr.choice(np.setdiff1d(np.arange(10, n - 10), np.concatenate([drop, drop - 1])), 3, replace=False)
    spur = opr_t[spur_i] + rng_opr.uniform(0.2, 0.8, 3) * (opr_t[spur_i + 1] - opr_t[spur_i])
    opr_t = np.sort(np.concatenate([np.delete(opr_t, drop), spur]))
    np.savetxt(os.path.join(data_dir, "opr.csv"), opr_t, fmt="%.9f", header="opr_time_s", comments="")

    def hermite(i, tq):
        k = np.clip(np.floor(tq / DT).astype(int), 0, nt - 2)
        s = (tq - t[k]) / DT
        h00 = 2*s**3 - 3*s**2 + 1; h10 = s**3 - 2*s**2 + s
        h01 = -2*s**3 + 3*s**2;    h11 = s**3 - s**2
        return h00 * X[k, i] + h10 * DT * V[k, i] + h01 * X[k + 1, i] + h11 * DT * V[k + 1, i]

    fres = EO * f0 / np.sqrt(EO**2 - S)
    t_res = (fres / EO * 60 - RPM0) / (RPM1 - RPM0) * TDUR
    gap = (np.quantile(t_res, 0.3), np.quantile(t_res, 0.65))  # dropout inside resonance region
    alpha = np.deg2rad(probe_nom + probe_err)
    for p in range(NPR):
        toas = []
        for i in range(B):
            m = np.arange(0, n_rev + 2)
            target = alpha[p] + 2 * np.pi * m + 2 * np.pi * i / B
            target = target[(target > Phi[0] + 0.1) & (target < Phi[-1] - 0.1)]
            # Phi may decrease? no: speed always positive, Phi monotonic
            tq = np.interp(target, Phi, t)
            for _ in range(4):
                Omq = np.interp(tq, t, Om)
                stat = stagger[i] + untwist[i] * (Omq / om_mid) ** 2
                defl = hermite(i, tq) + stat
                res = np.interp(tq, t, Phi) + defl / R - target
                tq = tq - res / (2 * np.pi * Omq)
            toas.append(tq)
        toa = np.sort(np.concatenate(toas)) + JIT_TOA * rng.standard_normal(sum(len(a) for a in toas))
        keep = rng.random(len(toa)) > 0.0007
        if p == gap_probe - 1:
            keep &= ~((toa > gap[0]) & (toa < gap[1]))
        toa = toa[keep]
        spur = rng.uniform(toa[0], toa[-1], int(0.0002 * len(toa)) + 3)
        toa = np.sort(np.concatenate([toa, spur]))
        np.savetxt(os.path.join(data_dir, f"probe_{p+1}.csv"), toa, fmt="%.9f", header="toa_s", comments="")

    rig = {
        "n_blades": B,
        "tip_radius_m": R,
        "southwell_coefficient": S,
        "probe_angles_deg": [round(float(a), 1) for a in probe_nom],
        "probe_angle_tolerance_deg": 0.1,
        "nominal_speed_rpm": {"start": RPM0, "end": RPM1},
        "nominal_duration_s": TDUR,
    }
    with open(os.path.join(data_dir, "rig.json"), "w") as f:
        json.dump(rig, f, indent=2)

    truth = {"engine_order": EO,
             "blades": {str(i + 1): {"resonance_frequency_Hz": float(fres[i]),
                                      "damping_ratio": float(zeta[i]),
                                      "steady_resonant_amplitude_mm": float(amp[i] / (2 * zeta[i]) * 1e3),
                                      "f0_Hz": float(f0[i])} for i in range(B)},
             "probe_error_deg": probe_err.tolist(), "seed": seed}
    os.makedirs(os.path.dirname(truth_path), exist_ok=True)
    with open(truth_path, "w") as f:
        json.dump(truth, f, indent=1)
    print(f"{name}: B={B} EO={EO} probes={probe_nom.tolist()} revs={n_rev} "
          f"res rpm {(fres/EO*60).min():.0f}-{(fres/EO*60).max():.0f} "
          f"peak mm {(amp/(2*zeta)).min()*1e3:.2f}-{(amp/(2*zeta)).max()*1e3:.2f}")


def main():
    only = os.environ.get("BTT_ONLY")
    for name in CONFIGS:
        if only and name not in only.split(","):
            continue
        data_dir = os.path.join(ROOT, "tests", "heldout", name)
        generate(name, data_dir, os.path.join(ROOT, "tests", "truth", f"{name}.json"))
        if name == "example":  # the shipped copy the agent sees
            import shutil
            dst = os.path.join(ROOT, "environment", "data", "example")
            shutil.rmtree(dst, ignore_errors=True)
            shutil.copytree(data_dir, dst)


if __name__ == "__main__":
    main()
