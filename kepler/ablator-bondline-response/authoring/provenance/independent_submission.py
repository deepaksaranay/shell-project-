#!/usr/bin/env python3
"""Structurally independent submission (written separately from the
answer-key generator and from solution/solve.py), scored through the real
verifier. Different at every stage:

  * grid: material-fixed (Lagrangian) cell-centred finite volumes on a
    uniform grid; no grid-motion convection terms. Recession consumes the
    first cell, which is merged into its neighbour (mass and energy
    conserving) when it becomes thin;
  * time integration: own fixed-step backward Euler with Picard iterations;
    decomposition updated implicitly cell by cell, energy solved as a
    linearised tridiagonal system in enthalpy form;
  * surface: T_w from the surface energy balance against the half-cell
    conduction resistance, linearised in the first cell temperature;
  * bondline: interface temperature from the series-resistance face value
    between ablator and aluminium cells.

Usage: python3 independent_submission.py <data_dir> <output_json>
"""
import json
import math
import os
import sys

import numpy as np
from scipy.linalg import solve_banded
from scipy.optimize import brentq

SIG = 5.670374419e-8
DX0 = 0.05e-3
DT = 0.05


def table(d, name):
    return np.loadtxt(os.path.join(d, name), delimiter=",", skiprows=1)


def enthalpy_fn(T_tab, cp_tab, hf):
    """h(T) = hf + int_298^T cp, cp linear between table points, held outside."""
    knots_H = np.concatenate([[0.0], np.cumsum(np.diff(T_tab) * (cp_tab[1:] + cp_tab[:-1]) / 2)])
    def H(T):
        T = np.atleast_1d(np.asarray(T, float))
        Tc = np.clip(T, T_tab[0], T_tab[-1])
        k = np.minimum(np.searchsorted(T_tab, Tc, side="right") - 1, len(T_tab) - 2)
        u = Tc - T_tab[k]
        g = (cp_tab[k + 1] - cp_tab[k]) / (T_tab[k + 1] - T_tab[k])
        val = knots_H[k] + cp_tab[k] * u + g * u * u / 2
        val += (T - Tc) * np.where(T < T_tab[0], cp_tab[0], cp_tab[-1])
        return val
    ref = H(298.0)[0]
    return lambda T: hf + H(T) - ref


def main(d, out_path):
    case = json.load(open(os.path.join(d, "case.json")))
    ab = case["ablator"]; dec = ab["decomposition"]; Gm = dec["Gamma"]; comp = dec["components"]
    pa, pb = comp["A"], comp["B"]; rhoC = comp["C"]["rho_kg_per_m3"]
    rho_v = Gm * (pa["rho_v_kg_per_m3"] + pb["rho_v_kg_per_m3"]) + (1 - Gm) * rhoC
    rho_c = Gm * (pa["rho_r_kg_per_m3"] + pb["rho_r_kg_per_m3"]) + (1 - Gm) * rhoC
    vt = table(d, ab["virgin"]["properties_file"]); ct = table(d, ab["char"]["properties_file"])
    h_v = enthalpy_fn(vt[:, 0], vt[:, 1], ab["virgin"]["formation_enthalpy_J_per_kg_at_298K"])
    h_c = enthalpy_fn(ct[:, 0], ct[:, 1], ab["char"]["formation_enthalpy_J_per_kg_at_298K"])
    gt = table(d, ab["pyrolysis_gas"]["enthalpy_file"])
    st = table(d, case["substrate"]["properties_file"]); rho_s = case["substrate"]["density_kg_per_m3"]
    h_s = enthalpy_fn(st[:, 0], st[:, 1], 0.0)
    th = table(d, case["surface"]["thermochemistry_file"])
    tT = np.unique(th[:, 0]); tG = np.unique(th[:, 1])
    tBc = th[:, 2].reshape(len(tT), len(tG)); tHw = th[:, 3].reshape(len(tT), len(tG))
    traj = table(d, "trajectory.csv")
    L = ab["thickness_mm"] * 1e-3; Ls = case["substrate"]["thickness_mm"] * 1e-3
    Tsink = case["surface"]["radiation_sink_temperature_K"]
    times = case["output_times_s"]

    def hg(T): return np.interp(T, gt[:, 0], gt[:, 1])
    def dhg(T):
        e = 0.5
        return (hg(T + e) - hg(T - e)) / (2 * e)

    def lookup(tab, T, g):
        T = min(max(T, tT[0]), tT[-1]); g = min(max(g, tG[0]), tG[-1])
        i = min(int(np.searchsorted(tT, T, side="right")) - 1, len(tT) - 2)
        j = min(int(np.searchsorted(tG, g, side="right")) - 1, len(tG) - 2)
        a = (T - tT[i]) / (tT[i + 1] - tT[i]); b = (g - tG[j]) / (tG[j + 1] - tG[j])
        return (1 - a) * (1 - b) * tab[i, j] + a * (1 - b) * tab[i + 1, j] + (1 - a) * b * tab[i, j + 1] + a * b * tab[i + 1, j + 1]

    # ---- cells: ablator (material-fixed), then substrate
    na = int(round(L / DX0)); ns = 12
    xa_faces = np.linspace(0.0, L, na + 1)
    Ws = Ls / ns
    W = np.diff(xa_faces).copy()                 # ablator cell widths (first one shrinks)
    MA = W * pa["rho_v_kg_per_m3"]; MB = W * pb["rho_v_kg_per_m3"]
    Ta = np.full(na, case["initial_temperature_K"]); Tsb = np.full(ns, case["initial_temperature_K"])
    first = 0                                    # index of first active ablator cell
    s = 0.0; mloss = 0.0

    def props(T, rho):
        y = np.clip(rho_v * (rho - rho_c) / (rho * (rho_v - rho_c)), 0, 1)
        cp = y * np.interp(T, vt[:, 0], vt[:, 1]) + (1 - y) * np.interp(T, ct[:, 0], ct[:, 1])
        k = y * np.interp(T, vt[:, 0], vt[:, 2]) + (1 - y) * np.interp(T, ct[:, 0], ct[:, 2])
        h = y * h_v(T) + (1 - y) * h_c(T)
        return y, cp, k, h

    def decompose(T, rA0, rB0):
        """Implicit Euler for each component density, vectorised Newton."""
        res = []
        for P, r0 in ((pa, rA0), (pb, rB0)):
            kk = P["k_per_s"] * np.exp(-P["E_over_R_K"] / T); rv, rr, ps = P["rho_v_kg_per_m3"], P["rho_r_kg_per_m3"], P["psi"]
            r = r0.copy()
            for _ in range(60):
                f = np.clip((r - rr) / rv, 0, None)
                F = r - r0 + DT * kk * rv * f ** ps
                dF = 1 + DT * kk * ps * f ** (ps - 1)
                step = F / dF
                r = np.maximum(r - step, rr)
                if np.max(np.abs(step)) < 1e-10:
                    break
            res.append(r)
        return res

    def surface_flux(t, Tw, mg, y_w, h_wall_solid):
        C0 = np.interp(t, traj[:, 0], traj[:, 1]); hr = np.interp(t, traj[:, 0], traj[:, 2]); qr = np.interp(t, traj[:, 0], traj[:, 3])
        def g(C):
            bg = mg / C; B = bg + lookup(tBc, Tw, bg)
            return C - C0 * (math.log(1 + B) / B if B > 1e-12 else 1.0)
        C = C0 if g(C0) <= 0 else brentq(g, 1e-14, C0, xtol=1e-16, rtol=1e-13)
        bg = mg / C
        mc = lookup(tBc, Tw, bg) * C; hw = lookup(tHw, Tw, bg)
        eps = y_w * np.interp(Tw, vt[:, 0], vt[:, 3]) + (1 - y_w) * np.interp(Tw, ct[:, 0], ct[:, 3])
        q = (C * (hr - hw) + mg * hg(Tw) + mc * h_wall_solid(Tw) - (mg + mc) * hw
             + eps * qr - eps * SIG * (Tw ** 4 - Tsink ** 4))
        return q, mc

    snaps = []; peak = 0.0
    nsteps = int(round(times[-1] / DT))
    out_idx = {int(round(tt / DT)): tt for tt in times}
    t = 0.0
    for n in range(1, nsteps + 1):
        t = n * DT
        act = slice(first, na)
        Wn = W[act].copy(); rA0 = MA[act] / Wn; rB0 = MB[act] / Wn
        rho0 = Gm * (rA0 + rB0) + (1 - Gm) * rhoC
        _, _, _, h0 = props(Ta[act], rho0)
        E0 = Wn * rho0 * h0
        Es0 = rho_s * Ws * h_s(Tsb)
        T = Ta[act].copy(); Tsub = Tsb.copy()
        sdot = 0.0
        for it in range(80):
            rA, rB = decompose(T, rA0, rB0)
            rho = Gm * (rA + rB) + (1 - Gm) * rhoC
            prod = -(rho - rho0) / DT * Wn          # gas generation per cell (kg/m2 s); W change below handled via consumption
            mf = np.concatenate([np.cumsum(prod[::-1])[::-1], [0.0]])  # mf[i]: through left face of cell i
            mg = mf[0]
            y, cp, k, h = props(T, rho)
            # surface: q(Tw) = K (Tw - T1), K = k1 / (W1/2)
            hsurf = lambda Tw, y1=y[0]: float(y1 * h_v(Tw)[0] + (1 - y1) * h_c(Tw)[0])
            def solve_w(T1, Wf):
                K = k[0] / (Wf / 2)
                fz = lambda Tw: surface_flux(t, Tw, mg, y[0], hsurf)[0] - K * (Tw - T1)
                Tw = brentq(fz, 150.0, 6000.0, xtol=1e-8)
                q, mc = surface_flux(t, Tw, mg, y[0], hsurf)
                return Tw, q, mc
            Tw, qk, mc = solve_w(T[0], Wn[0])
            sdot = mc / rho[0]
            W1 = Wn[0] - sdot * DT                     # surface cell after consumption
            Tw2, qk2, _ = solve_w(T[0] + 0.5, Wn[0])
            dq = (qk2 - qk) / 0.5                      # d(conducted flux)/dT1
            Wt = Wn.copy(); Wt[0] = W1
            # energy residual linearised in T: Wt*rho*h(T) - E0 + DT*(div fluxes) + consumption = 0
            m = len(T); N = m + ns
            kf = 0.5 * (k[1:] + k[:-1]); dxc = 0.5 * (Wt[1:] + Wt[:-1])
            Kf = kf / dxc
            # interface: series resistance ablator half cell + substrate half cell
            ksb = np.interp(Tsub, st[:, 0], st[:, 2])
            Ki = 1.0 / (Wt[-1] / 2 / k[-1] + Ws / 2 / ksb[0])
            Kss = 1.0 / (Ws / 2 / ksb[1:] + Ws / 2 / ksb[:-1])
            cps = np.interp(Tsub, st[:, 0], st[:, 1])
            cap = np.concatenate([Wt * rho * cp, rho_s * Ws * cps])
            Ecur = np.concatenate([Wt * rho * h, rho_s * Ws * h_s(Tsub)])
            Eold = np.concatenate([E0, Es0])
            Tall = np.concatenate([T, Tsub])
            # residual R = (E - Eold)/DT - net inflow ; Newton-like linearisation
            R = (Ecur - Eold) / DT
            J = np.zeros((3, N)); J[1] = cap / DT
            # conduction faces
            Kall = np.concatenate([Kf, [Ki], Kss])
            fl = Kall * (Tall[1:] - Tall[:-1])          # heat flow from i+1 into i
            R[:-1] -= fl; R[1:] += fl
            J[1, :-1] += Kall; J[0, 1:] -= Kall          # dR_i/dT_{i+1} stored at upper diag
            J[1, 1:] += Kall; J[2, :-1] -= Kall
            # pyrolysis gas: enters cell i from i+1 with h_g(T_{i+1}), leaves to i-1 with h_g(T_i)
            gin = mf[1:m] * hg(T[1:])                    # into cell i (i<m-1) from below
            R[:m - 1] -= gin
            J[0, 1:m] -= mf[1:m] * dhg(T[1:])
            gout = mf[1:m] * hg(T[1:])                   # same stream leaves cell i+1
            R[1:m] += gout
            J[1, 1:m] += mf[1:m] * dhg(T[1:])
            # first cell: gas leaves at the wall temperature through the surface, conducted flux, consumption
            # net into the first cell: conducted flux q_k (the SEB value, which contains
            # mdot_g*h_g(T_w) + mdot_c*h_s(T_w)) minus the gas leaving at T_w minus
            # the consumed surface slice, which leaves with h_s(T_w)
            R[0] += -qk + mg * hg(Tw) + mc * hsurf(Tw)
            J[1, 0] += -dq
            dT = solve_banded((1, 1), J, -R)
            Tall = Tall + dT
            T = Tall[:m]; Tsub = Tall[m:]
            if np.max(np.abs(dT)) < 1e-7:
                break
        # accept step
        rA, rB = decompose(T, rA0, rB0)
        Wn[0] = W1
        MA[act] = rA * Wn; MB[act] = rB * Wn; W[act] = Wn
        Ta[act] = T; Tsb = Tsub
        s += sdot * DT; mloss += (mg + mc) * DT
        # merge thin surface cell into its neighbour (conserves mass and energy)
        if W[first] < 0.4 * DX0:
            i, j = first, first + 1
            y1, _, _, h1 = props(Ta[i:i + 2], (Gm * (MA[i:i + 2] + MB[i:i + 2]) / W[i:i + 2] + (1 - Gm) * rhoC))
            Etot = (Gm * (MA[i] + MB[i]) + (1 - Gm) * rhoC * W[i]) * h1[0] + (Gm * (MA[j] + MB[j]) + (1 - Gm) * rhoC * W[j]) * h1[1]
            MA[j] += MA[i]; MB[j] += MB[i]; W[j] += W[i]
            rj = Gm * (MA[j] + MB[j]) / W[j] + (1 - Gm) * rhoC
            Mj = rj * W[j]
            Tj = Ta[j]
            for _ in range(50):
                yj, cpj, _, hj = props(np.array([Tj]), np.array([rj]))
                f = Mj * hj[0] - Etot
                Tj -= f / (Mj * cpj[0])
                if abs(f) < 1e-6:
                    break
            Ta[j] = Tj
            first += 1
        # bondline: face temperature between last ablator cell and first substrate cell
        rl = Gm * (MA[-1] + MB[-1]) / W[-1] + (1 - Gm) * rhoC
        kl = props(Ta[-1:], np.array([rl]))[2][0]; ks0 = np.interp(Tsb[0], st[:, 0], st[:, 2])
        ga, gs = kl / (W[-1] / 2), ks0 / (Ws / 2)
        Tb = (ga * Ta[-1] + gs * Tsb[0]) / (ga + gs)
        peak = max(peak, Tb)
        if n in out_idx:
            act = slice(first, na)
            rho = Gm * (MA[act] + MB[act]) / W[act] + (1 - Gm) * rhoC
            tau = (rho - rho_c) / (rho_v - rho_c)
            xc = np.cumsum(W[act]) - W[act] / 2       # cell centres below the current surface
            def front(lv):
                kk = np.nonzero(tau >= lv)[0]
                if len(kk) == 0:
                    return xc[-1]
                kk = kk[0]
                if kk == 0:
                    return 0.0
                return xc[kk - 1] + (lv - tau[kk - 1]) / (tau[kk] - tau[kk - 1]) * (xc[kk] - xc[kk - 1])
            snaps.append({"t_s": out_idx[n], "recession_mm": s * 1e3, "surface_temperature_K": float(Tw),
                          "bondline_temperature_K": float(Tb), "pyrolysis_front_depth_mm": float(front(0.98) * 1e3),
                          "char_front_depth_mm": float(front(0.02) * 1e3), "mass_loss_kg_per_m2": float(mloss)})
            print(snaps[-1], flush=True)
    out = {"snapshots": snaps, "peak_bondline_temperature_K": float(peak)}
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
