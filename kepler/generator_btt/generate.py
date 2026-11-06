#!/usr/bin/env python3
"""Generate shipped data + sealed truth for the blade-tip-timing task.

NOT part of the bundle. Run: python3 kepler/generator_btt/generate.py

Each blade i is an uncoupled SDOF oscillator (tip, circumferential):
    x'' + 2 z_i wn_i(t) x' + wn_i(t)^2 x = wn_i(t)^2 * a_i * cos(EO*Phi(t) - EO*2pi(i-1)/B + psi)
    fn_i(Omega)^2 = f0_i^2 + S * Omega^2      (Southwell, Omega in rev/s)
integrated through a slow run-up. Tip arrival times at casing probes are
found from theta_i(t) + x_i(t)/R = alpha_p.
"""
import json
import os

import numpy as np

ROOT = os.path.join(os.path.dirname(__file__), "..", "blade-tip-timing")
DATA = os.path.join(ROOT, "environment", "data")
TESTS = os.path.join(ROOT, "tests")
rng = np.random.default_rng(7710)

B = 24
R = 0.2475               # tip radius, m
EO = 6
S = 4.2                  # Southwell coefficient (stated)
PROBE_NOM = np.array([0.0, 14.0, 48.5, 56.5, 153.0])  # deg, stated nominal
PROBE_ERR = rng.uniform(-0.08, 0.08, 5)         # unknown install error, deg
RPM0, RPM1, TDUR = 6550.0, 7550.0, 44.0
DT = 1.0 / 25000.0
JIT_TOA = 0.30e-6
JIT_OPR = 0.15e-6

f0 = 662.0 * (1 + 0.012 * rng.standard_normal(B))
zeta = rng.uniform(0.0025, 0.0060, B)
amp = 3.2e-6 * (1 + 0.15 * rng.standard_normal(B))   # static-equivalent forcing, m
stagger = 0.4e-3 * rng.standard_normal(B)            # static tip position error, m
untwist = 0.25e-3 * (1 + 0.3 * rng.standard_normal(B))  # static drift at 7000 rpm, m
psi = rng.uniform(0, 2 * np.pi)


def speed_hz(t):
    lin = (RPM0 + (RPM1 - RPM0) * t / TDUR) / 60.0
    return lin + 0.04 * np.sin(2 * np.pi * 0.37 * t) + 0.02 * np.sin(2 * np.pi * 1.9 * t + 1.0)


def main():
    os.makedirs(DATA, exist_ok=True)
    nt = int(TDUR / DT) + 1
    t = np.arange(nt) * DT
    Om = speed_hz(t)
    # rotor angle (rad), trapezoid
    Phi = np.concatenate([[0.0], np.cumsum(0.5 * (Om[1:] + Om[:-1]) * DT)]) * 2 * np.pi
    Phi += 0.3

    # integrate SDOFs with RK4, coefficients evaluated at sub-steps by interpolation
    blade_ph = EO * 2 * np.pi * np.arange(B) / B
    w0sq = (2 * np.pi * f0) ** 2
    def rhs(Phi_v, Om_v, x, v):
        wn2 = w0sq + S * (2 * np.pi * Om_v) ** 2
        wn = np.sqrt(wn2)
        f = wn2 * amp * np.cos(EO * Phi_v - blade_ph + psi)
        return v, f - 2 * zeta * wn * v - wn2 * x
    X = np.zeros((nt, B), np.float64)
    V = np.zeros((nt, B), np.float64)
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

    # OPR: Phi = 2 pi n
    n_rev = int(np.floor((Phi[-1] - 0.05) / (2 * np.pi)))
    opr_t = np.interp(2 * np.pi * np.arange(1, n_rev + 1), Phi, t)
    opr_t += JIT_OPR * rng.standard_normal(len(opr_t))
    np.savetxt(os.path.join(DATA, "opr.csv"), opr_t, fmt="%.9f", header="opr_time_s", comments="")

    def hermite(i, tq):
        k = np.clip(np.floor(tq / DT).astype(int), 0, nt - 2)
        s = (tq - t[k]) / DT
        h00 = 2*s**3 - 3*s**2 + 1; h10 = s**3 - 2*s**2 + s
        h01 = -2*s**3 + 3*s**2;    h11 = s**3 - s**2
        return (h00 * X[k, i] + h10 * DT * V[k, i] + h01 * X[k + 1, i] + h11 * DT * V[k + 1, i])

    alpha = np.deg2rad(PROBE_NOM + PROBE_ERR)
    for p in range(len(alpha)):
        toas = []
        for i in range(B):
            ang0 = -2 * np.pi * i / B     # blade i trails blade 1 by i pitches
            # undeflected arrival: Phi + ang0 + static/R = alpha + 2 pi m
            m = np.arange(0, n_rev + 2)
            target = alpha[p] + 2 * np.pi * m - ang0
            ok = (target > Phi[0] + 0.1) & (target < Phi[-1] - 0.1)
            target = target[ok]
            tq = np.interp(target, Phi, t)
            for _ in range(4):
                Omq = np.interp(tq, t, Om)
                stat = stagger[i] + untwist[i] * (Omq / (7000 / 60)) ** 2
                defl = hermite(i, tq) + stat
                # condition: Phi(t) + ang0 + defl/R = alpha + 2 pi m
                Phiq = np.interp(tq, t, Phi)
                res = Phiq + defl / R - target
                tq = tq - res / (2 * np.pi * Omq)
            toas.append(tq)
        toa = np.sort(np.concatenate(toas))
        toa = toa + JIT_TOA * rng.standard_normal(len(toa))
        # dropouts and spurious triggers
        keep = rng.random(len(toa)) > 0.0007
        toa = toa[keep]
        spur = rng.uniform(toa[0], toa[-1], int(0.0002 * len(toa)) + 3)
        toa = np.sort(np.concatenate([toa, spur]))
        np.savetxt(os.path.join(DATA, f"probe_{p+1}.csv"), toa, fmt="%.9f",
                   header="toa_s", comments="")

    fres = EO * f0 / np.sqrt(EO**2 - S)
    truth = {"engine_order": EO,
             "blades": {str(i + 1): {"resonance_frequency_Hz": float(fres[i]),
                                      "damping_ratio": float(zeta[i]),
                                      "f0_Hz": float(f0[i])} for i in range(B)},
             "probe_error_deg": PROBE_ERR.tolist()}
    with open(os.path.join(TESTS, "expected.json"), "w") as f:
        json.dump(truth, f, indent=1)
    print("revs", n_rev, "fres range", fres.min(), fres.max(),
          "res rpm", (fres / EO * 60).min(), (fres / EO * 60).max())
    peak = amp / (2 * zeta)
    print("peak tip amp mm", peak.min() * 1e3, peak.max() * 1e3)


if __name__ == "__main__":
    main()
