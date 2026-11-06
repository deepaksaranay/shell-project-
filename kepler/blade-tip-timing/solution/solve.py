#!/usr/bin/env python3
"""Reference blade-tip-timing analysis.

1. Shaft angle Phi(t) from a cubic spline through the once-per-rev times.
2. Each probe arrival -> shaft angle -> blade assignment by nearest pitch slot;
   spurious triggers (far from any slot / duplicates) are discarded.
3. Tip deflection d = R * (slot angle - arrival shaft angle) (tip ahead of
   its slot arrives early).
4. Per revolution and blade, least-squares fit across the probes:
       d_p = c + a cos(EO alpha_p) + b sin(EO alpha_p)
   -> complex synchronous response Q = a - i b. EO is chosen as the integer
   whose fit leaves the smallest residual during the resonance.
5. Per blade, fit Q(Omega) with a centrifugally stiffened SDOF:
       Q = G wn^2 / (wn^2 - w^2 + 2 i zeta wn w) + c,  w = 2 pi EO Omega,
       wn^2 = (2 pi f0)^2 + S (2 pi Omega)^2
   (complex constant c absorbs probe-angle installation errors).
   f_res = EO * f0 / sqrt(EO^2 - S).
   This steady-state fit only seeds step 6: the run-up is fast compared with
   the blades' decay time, so the measured curve is a transient (peak shifted
   and broadened, ringing after the crossing). Fitting it as steady state
   biases both frequency and damping.
6. Per blade, fit a transient model to every probe sample: the slow complex
   envelope of the SDOF response, integrated under the measured speed
   history (Omega(t) from the OPR spline), plus per-probe constants and a
   quadratic common-mode drift (untwist). (f0, zeta) nonlinear; the complex
   gain and offsets linear (variable projection); 5-sigma outlier trimming.
"""
import argparse
import json
import os

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import least_squares

B = 24
R = 0.2475
S = 4.2
ALPHA = np.deg2rad([0.0, 14.0, 48.5, 56.5, 153.0])
NP = len(ALPHA)


def load(data):
    opr = np.loadtxt(os.path.join(data, "opr.csv"), skiprows=1)
    toas = [np.loadtxt(os.path.join(data, f"probe_{p}.csv"), skiprows=1) for p in range(1, NP + 1)]
    return opr, toas


def deflections(opr, toas):
    """Return D[rev, blade, probe] deflection (nan where missing), Omega[rev]
    (rev/s), arrival times Tm[rev, blade, probe] and the shaft-angle spline."""
    sp = CubicSpline(opr, 2 * np.pi * np.arange(len(opr)))
    nrev = len(opr) - 1
    D = np.full((nrev, B, NP), np.nan)
    Tm = np.full((nrev, B, NP), np.nan)
    pitch = 2 * np.pi / B
    for p, toa in enumerate(toas):
        toa = toa[(toa > opr[0]) & (toa < opr[-1])]
        ph = sp(toa)
        rev = np.floor(ph / (2 * np.pi)).astype(int)
        beta = ph - 2 * np.pi * rev
        rel = (beta - ALPHA[p]) / pitch
        slot = np.round(rel).astype(int)
        err = rel - slot
        # slots past 2*pi wrap into the next OPR revolution; the blade is the
        # slot modulo B and the response is quasi-steady across one rev
        blade = np.mod(slot, B)
        rev_b = rev
        d = -R * err * pitch
        good = np.abs(err) < 0.06
        seen = np.zeros((nrev, B), bool)
        for r, b_, dd, tt in zip(rev_b[good], blade[good], d[good], toa[good]):
            if 0 <= r < nrev:
                if seen[r, b_]:
                    D[r, b_, p] = np.nan  # duplicate trigger: ambiguous, drop
                else:
                    D[r, b_, p] = dd
                    Tm[r, b_, p] = tt
                    seen[r, b_] = True
    Om = 1.0 / np.diff(opr)
    return D, Om, Tm, sp


def per_rev_fit(D, eo):
    """Return complex synchronous response Q[rev, blade] and residual SS."""
    nrev = D.shape[0]
    Q = np.full((nrev, B), np.nan, complex)
    res = np.full((nrev, B), np.nan)
    A_full = np.column_stack([np.ones(NP), np.cos(eo * ALPHA), np.sin(eo * ALPHA)])
    finite = np.isfinite(D)
    full = finite.all(axis=2)
    Y = D[full]                                   # (n, NP)
    coef = Y @ np.linalg.pinv(A_full).T           # (n, 3)
    Q[full] = coef[:, 1] - 1j * coef[:, 2]
    res[full] = np.sum((Y - coef @ A_full.T) ** 2, axis=1)
    for r, b in zip(*np.where(~full & (finite.sum(axis=2) >= 3))):
        m = finite[r, b]
        c, *_ = np.linalg.lstsq(A_full[m], D[r, b, m], rcond=None)
        Q[r, b] = c[1] - 1j * c[2]
        res[r, b] = np.sum((D[r, b, m] - A_full[m] @ c) ** 2) * NP / m.sum()
    return Q, res


def zero_static(D, Om):
    """Remove a constant per (blade, probe). Speed-dependent untwist is common
    to all probes of a blade and drops out in the per-revolution constant; a
    probe-angle error is a constant per probe and maps to a constant complex
    offset, which the SDOF fit carries explicitly."""
    return D - np.nanmedian(D, axis=0, keepdims=True)


def clean_Q(Q, k=6.0, half=6):
    """Drop revolutions whose response jumps away from its neighbours
    (spurious trigger that survived slot gating)."""
    Q = Q.copy()
    n = Q.shape[0]
    for b in range(Q.shape[1]):
        q = Q[:, b]
        med = np.full(n, np.nan, complex)
        for r in range(n):
            seg = q[max(0, r - half):r + half + 1]
            seg = seg[np.isfinite(seg)]
            if len(seg) >= 3:
                med[r] = np.median(seg.real) + 1j * np.median(seg.imag)
        dev = np.abs(q - med)
        sig = 1.4826 * np.nanmedian(dev)
        q[dev > k * sig] = np.nan
    return Q


def choose_eo(D):
    """EO whose per-revolution fit explains the probe-to-probe scatter best."""
    total = np.nansum(D ** 2, axis=2)
    best = None
    for eo in range(1, 25):
        Q, res = per_rev_fit(D, eo)
        # weight revolutions by how much they vibrate (resonance region)
        wgt = np.nan_to_num(total)
        score = np.nansum(res * wgt) / np.nansum(total * wgt)
        if best is None or score < best[0]:
            best = (score, eo, Q)
    return best[1], best[2]


def fit_blade(Om, q, eo, use_southwell=True, half_window=2.5):
    """Variable-projection SDOF fit: (f0, zeta) nonlinear, (G, c) linear."""
    s = S if use_southwell else 0.0
    m = np.isfinite(q)
    Om, q = Om[m], q[m]
    sm = np.convolve(np.abs(q), np.ones(15) / 15, "same")
    Om_pk = Om[np.argmax(sm)]
    win = np.abs(Om - Om_pk) < half_window
    Om, q = Om[win], q[win]
    w = 2 * np.pi * eo * Om

    def resid(p):
        f0, lz = p
        z = np.exp(lz)
        wn2 = (2 * np.pi * f0) ** 2 + s * (2 * np.pi * Om) ** 2
        h = wn2 / (wn2 - w**2 + 2j * z * np.sqrt(wn2) * w)
        A = np.column_stack([h, np.ones_like(h)])
        coef, *_ = np.linalg.lstsq(A, q, rcond=None)
        r = A @ coef - q
        return np.concatenate([r.real, r.imag])

    f0c = Om_pk * np.sqrt(eo**2 - s)
    best = None
    for f0 in f0c + np.linspace(-2.0, 2.0, 21):
        for lz in np.log(np.geomspace(1e-3, 2e-2, 13)):
            c = np.sum(resid([f0, lz]) ** 2)
            if best is None or c < best[0]:
                best = (c, f0, lz)
    sol = least_squares(resid, [best[1], best[2]], x_scale=[0.1, 0.05])
    f0, z = sol.x[0], float(np.exp(sol.x[1]))
    return eo * f0 / np.sqrt(eo**2 - s), z


def envelope(tg, Om_g, f0, z, eo, s=S):
    """Slow complex envelope u(t) of the SDOF response to unit engine-order
    forcing, x = Re(u e^{i EO Phi}), under the measured speed history.
    u' = [wn^2 - (wn^2 - w^2 + 2 i z wn w) u] / (2 i w + 2 z wn), exact for
    piecewise-constant coefficients; starts from the quasi-static value."""
    w = 2 * np.pi * eo * Om_g
    wn2 = (2 * np.pi * f0) ** 2 + s * (2 * np.pi * Om_g) ** 2
    wn = np.sqrt(wn2)
    den = 2j * w + 2 * z * wn
    a = -(wn2 - w**2 + 2j * z * wn * w) / den
    b = wn2 / den
    h = np.diff(tg)
    ah = a[:-1] * h
    alpha = np.exp(ah)
    beta = b[:-1] * (alpha - 1) / a[:-1]
    logP = np.concatenate([[0.0], np.cumsum(ah)])
    P = np.exp(logP)
    u0 = wn2[0] / (wn2[0] - w[0] ** 2 + 2j * z * wn[0] * w[0])
    acc = np.concatenate([[0.0], np.cumsum(beta * np.exp(-logP[1:]))])
    return P * (u0 + acc)


def fit_blade_transient(D, Tm, Qc, sp, eo, b, f0_init, z_init, s=S, half=1.6):
    """Fit the transient SDOF model directly to every probe sample of blade b."""
    q = Qc[:, b]
    tr = np.nanmean(Tm[:, b, :], axis=1)
    sm = np.convolve(np.nan_to_num(np.abs(q)), np.ones(9) / 9, "same")
    t_pk = tr[np.nanargmax(sm)]
    rows = np.isfinite(q) & (np.abs(tr - t_pk) < half)
    d = D[rows, b, :]
    tt = Tm[rows, b, :]
    m = np.isfinite(d)
    y, ts = d[m], tt[m]
    pidx = np.nonzero(m)[1]
    tg = np.arange(ts.min() - 0.6, ts.max() + 0.01, 1e-3)
    Om_g = sp(tg, 1) / (2 * np.pi)
    Om_g = np.convolve(np.pad(Om_g, 25, mode="edge"), np.ones(51) / 51, "valid")
    rot = np.exp(1j * eo * ALPHA[pidx])
    tau = ts - t_pk
    base = [np.eye(NP)[pidx], tau[:, None], tau[:, None] ** 2]
    keep = np.ones(len(y), bool)

    def design(p):
        f0, lz = p
        u = envelope(tg, Om_g, f0, np.exp(lz), eo, s)
        ui = np.interp(ts, tg, u.real) + 1j * np.interp(ts, tg, u.imag)
        v = ui * rot
        return np.column_stack([v.real, -v.imag] + base)

    def resid(p):
        A = design(p)
        coef, *_ = np.linalg.lstsq(A[keep], y[keep], rcond=None)
        return (A @ coef - y)[keep]

    p0 = [f0_init, np.log(z_init)]
    for _ in range(2):
        sol = least_squares(resid, p0, x_scale=[0.05, 0.02])
        p0 = sol.x
        A = design(p0)
        coef, *_ = np.linalg.lstsq(A[keep], y[keep], rcond=None)
        r = A @ coef - y
        sig = 1.4826 * np.median(np.abs(r[keep]))
        keep = np.abs(r) < 5 * sig
    sol = least_squares(resid, p0, x_scale=[0.05, 0.02])
    f0, z = sol.x[0], float(np.exp(sol.x[1]))
    return eo * f0 / np.sqrt(eo**2 - s), z


def analyse(data, use_southwell=True, method="transient"):
    opr, toas = load(data)
    D, Om, Tm, sp = deflections(opr, toas)
    D = zero_static(D, Om)
    eo, Q = choose_eo(D)
    Q = clean_Q(Q)
    s = S if use_southwell else 0.0
    out = {"engine_order": int(eo), "blades": {}}
    for b in range(B):
        fres, z = fit_blade(Om, Q[:, b], eo, use_southwell)
        if method == "transient":
            f0 = fres * np.sqrt(eo**2 - s) / eo
            fres, z = fit_blade_transient(D, Tm, Q, sp, eo, b, f0, z, s)
        out["blades"][str(b + 1)] = {"resonance_frequency_Hz": float(fres), "damping_ratio": float(z)}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="/app/data")
    ap.add_argument("--out", default="/app/output/blade_modal_results.json")
    a = ap.parse_args()
    res = analyse(a.data)
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w") as f:
        json.dump(res, f, indent=2)


if __name__ == "__main__":
    main()
