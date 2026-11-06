#!/usr/bin/env python3
"""Reference identification of fluid-film bearing K, C, M from shaker records.

Steps
1. Convert voltages to engineering units; rotate the 45/135 deg probe pair
   into x/y housing-relative-to-shaft displacement (probe output is -gap).
2. Extract complex amplitudes at the multisine tones with a Hann window
   (tones are bin-centred, so the window only suppresses leakage from the
   unsynchronised 1X/2X runout lines; tones near 1X/2X are dropped anyway).
3. Dry baseline: (F - Ms*a) = Hs(w) * x_abs, x_abs = -a/w^2, two tests give
   the 2x2 support impedance at each tone.
4. Lubricated, rotating: (F - Ms*a - Hs*x_abs) = Hb(w) * d_rel, where d_rel
   comes from the probes. The support acts on ABSOLUTE housing motion, the
   film on motion RELATIVE to the (moving) shaft.
5. Least-squares fit Re Hb = K - w^2 M, Im Hb = w C per element.
"""
import argparse
import csv
import json
import os

import numpy as np

LC_SENS = 1.124e-3
ACC_SENS = 0.100 / 9.80665
PROX_SENS = 7.874e3
MS = 18.5
TONES = np.arange(20.0, 321.0, 10.0)
ANG = np.deg2rad([45.0, 135.0])


def load(path):
    a = np.loadtxt(path, delimiter=",", skiprows=1)
    return {
        "t": a[:, 0], "fx": a[:, 1] / LC_SENS, "fy": a[:, 2] / LC_SENS,
        "ax": a[:, 3] / ACC_SENS, "ay": a[:, 4] / ACC_SENS,
        "p1": a[:, 5], "p2": a[:, 6], "kp": a[:, 7],
    }


def shaft_speed_hz(r):
    t, kp = r["t"], r["kp"]
    thr = 0.5 * (kp.max() + kp.min())
    low = kp < thr
    edges = np.where(low[1:] & ~low[:-1])[0] + 1
    if len(edges) < 3:
        return None
    # regression of edge time on revolution index averages out sample jitter
    period = np.polyfit(np.arange(len(edges)), t[edges], 1)[0]
    return 1.0 / period


def probes_to_xy(r, opts):
    g1 = -(r["p1"] - r["p1"].mean()) / PROX_SENS
    g2 = -(r["p2"] - r["p2"].mean()) / PROX_SENS
    if opts.get("no_probe_sign"):
        g1, g2 = -g1, -g2
    if opts.get("no_rotation"):
        return g1, g2
    R = np.array([[np.cos(ANG[0]), np.sin(ANG[0])], [np.cos(ANG[1]), np.sin(ANG[1])]])
    xy = np.linalg.solve(R, np.vstack([g1, g2]))
    return xy[0], xy[1]


def tone_amps(x, fs, opts):
    n = len(x)
    w = np.ones(n) if opts.get("rect") else np.hanning(n)
    X = np.fft.rfft((x - x.mean()) * w)
    idx = np.round(TONES * n / fs).astype(int)
    return X[idx]


def spectra(r, opts):
    fs = 1.0 / (r["t"][1] - r["t"][0])
    dx, dy = probes_to_xy(r, opts)
    out = {}
    for k, sig in [("fx", r["fx"]), ("fy", r["fy"]), ("ax", r["ax"]), ("ay", r["ay"]),
                   ("dx", dx), ("dy", dy)]:
        out[k] = tone_amps(sig, fs, opts)
    return out


def mat(sx, sy, a, b):
    """2x2 per tone: rows = channel (x,y), cols = test (x-shaken, y-shaken)."""
    return np.stack([np.stack([sx[a], sy[a]], -1), np.stack([sx[b], sy[b]], -1)], -2)


def fit_kcm(H, w, keep):
    K = np.zeros((2, 2)); C = np.zeros((2, 2)); M = np.zeros((2, 2))
    ww = w[keep]
    Ar = np.column_stack([np.ones_like(ww), -ww**2])
    for i in range(2):
        for j in range(2):
            h = H[keep, i, j]
            (K[i, j], M[i, j]), *_ = np.linalg.lstsq(Ar, h.real, rcond=None)
            C[i, j] = np.sum(ww * h.imag) / np.sum(ww**2)
    return K, C, M


def identify(data_dir, opts=None):
    opts = opts or {}
    w = 2 * np.pi * TONES
    bx = spectra(load(os.path.join(data_dir, "baseline_dry_x.csv")), opts)
    by = spectra(load(os.path.join(data_dir, "baseline_dry_y.csv")), opts)
    F0, A0 = mat(bx, by, "fx", "fy"), mat(bx, by, "ax", "ay")
    X0 = -A0 / (w**2)[:, None, None]
    Hs = np.linalg.solve(np.swapaxes(X0, 1, 2), np.swapaxes(F0 - MS * A0, 1, 2))
    Hs = np.swapaxes(Hs, 1, 2)  # (F - Ms A) X^-1
    if opts.get("no_baseline"):
        Hs = np.zeros_like(Hs)

    results = {}
    with open(os.path.join(data_dir, "test_log.csv")) as f:
        rows = list(csv.DictReader(f))
    for row in rows:
        rx = load(os.path.join(data_dir, row["x_shaker_file"]))
        ry = load(os.path.join(data_dir, row["y_shaker_file"]))
        fr = shaft_speed_hz(rx) or float(row["nominal_speed_rpm"]) / 60.0
        sx, sy = spectra(rx, opts), spectra(ry, opts)
        F, A, D = mat(sx, sy, "fx", "fy"), mat(sx, sy, "ax", "ay"), mat(sx, sy, "dx", "dy")
        Xabs = -A / (w**2)[:, None, None]
        if opts.get("relative_for_support"):
            Xabs = D
        if opts.get("diag_force"):
            F = F * np.eye(2)[None]
        rhs = F - (0 if opts.get("no_inertia") else MS) * A - Hs @ Xabs
        Hb = np.swapaxes(np.linalg.solve(np.swapaxes(D, 1, 2), np.swapaxes(rhs, 1, 2)), 1, 2)
        keep = np.ones(len(TONES), bool)
        if not opts.get("keep_all"):
            for h in (1, 2, 3):
                keep &= np.abs(TONES - h * fr) > 2.0
        K, C, M = fit_kcm(Hb, w, keep)
        if opts.get("no_mass"):
            K = np.array([[np.mean(Hb[keep, i, j].real) for j in range(2)] for i in range(2)])
        results[row["condition"]] = {"K": K.tolist(), "C": C.tolist(), "M": M.tolist()}
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="/app/data")
    ap.add_argument("--out", default="/app/output/bearing_coefficients.json")
    a = ap.parse_args()
    res = identify(a.data)
    out = {c: {"K_N_per_m": v["K"], "C_Ns_per_m": v["C"]} for c, v in res.items()}
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w") as f:
        json.dump(out, f, indent=2)


if __name__ == "__main__":
    main()
