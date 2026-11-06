"""Method-of-characteristics water-hammer model of the rising main in
system.json, used to evaluate a candidate air-vessel design.

Fixed time step dt; each pipe gets an integer number of reaches with its wave
speed adjusted (|adjustment| reported); quasi-steady Darcy friction; pump
rundown from the homologous head/torque polynomials with trapezoidal speed
integration; ideal check valve; air vessel with polytropic gas, liquid level,
and differential connection losses. Pressures are checked at every
computational section (profile linearly interpolated).
"""
import json
import math

import numpy as np


def load_system(path):
    with open(path) as f:
        return json.load(f)


def simulate(sysd, design, dt=0.005, record=False):
    fl = sysd["fluid"]; g = fl["gravity_m_s2"]; Hatm = fl["atmospheric_pressure_head_m"]
    pump = sysd["pump"]; av = sysd["air_vessel"]
    Hs = sysd["suction_reservoir"]["water_level_m"]
    Hres = sysd["discharge_reservoir"]["water_level_m"]
    QR, HR = pump["rated_flow_m3_s"], pump["rated_head_m"]
    TR, I = pump["rated_torque_N_m"], pump["polar_moment_of_inertia_kg_m2"]
    wR = pump["rated_speed_rpm"] * 2 * math.pi / 60
    h0, h1, h2 = (pump["head_coefficients"][k] for k in ("h0", "h1", "h2"))
    b0, b1, b2 = (pump["torque_coefficients"][k] for k in ("b0", "b1", "b2"))
    chain = np.array(sysd["profile"]["chainage_m"]); elev = np.array(sysd["profile"]["elevation_m"])

    Vt = float(design["vessel_total_volume_m3"])
    Vg0 = float(design["initial_gas_volume_m3"])
    din = float(design["inflow_orifice_diameter_m"])
    n = av["gas_polytropic_exponent"]
    Dv = (4 * Vt / (av["height_to_diameter_ratio"] * math.pi)) ** (1 / 3)
    Av = math.pi * Dv**2 / 4
    zvb = av["bottom_elevation_m"]
    Ac = math.pi * av["outflow_connection"]["diameter_m"] ** 2 / 4
    Cout = av["outflow_connection"]["loss_coefficient"] / (2 * g * Ac**2)
    Ain = math.pi * din**2 / 4
    Cin = av["inflow_orifice"]["loss_coefficient"] / (2 * g * Ain**2)

    def pump_h(Qp, a):
        v = Qp / QR
        return HR * (h0 * a * a + h1 * a * v + h2 * v * v)

    def pump_T(Qp, a):
        v = Qp / QR
        return TR * (b0 * a * a + b1 * a * v + b2 * v * v)

    pipes = sysd["pipes"]
    A = [math.pi * p["diameter_m"] ** 2 / 4 for p in pipes]
    # steady state
    def system_resid(Q):
        hf = sum(p["darcy_friction_factor"] * p["length_m"] / p["diameter_m"] * (Q / Ai) ** 2 / (2 * g)
                 for p, Ai in zip(pipes, A))
        return Hs + pump_h(Q, 1.0) - Hres - hf
    lo, hi = 1e-6, 3 * QR
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if system_resid(mid) > 0:
            lo = mid
        else:
            hi = mid
    Q0 = 0.5 * (lo + hi)

    # grids
    N, B, R, x0, adj = [], [], [], [], []
    start = 0.0
    for p, Ai in zip(pipes, A):
        Ni = max(1, int(round(p["length_m"] / (p["wave_speed_m_s"] * dt))))
        dx = p["length_m"] / Ni
        a_adj = dx / dt
        adj.append(a_adj / p["wave_speed_m_s"] - 1)
        N.append(Ni)
        B.append(a_adj / (g * Ai))
        R.append(p["darcy_friction_factor"] * dx / (2 * g * p["diameter_m"] * Ai**2))
        x0.append(start + np.arange(Ni + 1) * dx)
        start += p["length_m"]
    z = [np.interp(x, chain, elev) for x in x0]

    H, Q = [], []
    Hcur = Hs + pump_h(Q0, 1.0)
    for i in range(3):
        Hi = Hcur - R[i] * Q0**2 * np.arange(N[i] + 1)
        H.append(Hi); Q.append(np.full(N[i] + 1, Q0))
        Hcur = Hi[-1]

    # vessel initial state
    zs = zvb + (Vt - Vg0) / Av
    HP0 = H[0][0]
    Hg = HP0 + Hatm - zs
    Cgas = Hg * Vg0**n
    Vg, Qv, Qp, a = Vg0, 0.0, Q0, 1.0
    cv_closed = False

    min_abs = min(float(np.min(H[i] - z[i] + Hatm)) for i in range(3))
    max_gauge = max(float(np.max(H[i] - z[i])) for i in range(3))
    max_gas = Vg / Vt
    where_min = where_max = (0.0, 0.0)
    rec = {"t": [], "Qp": [], "Qv": [], "HP": [], "alpha": [], "Vg": []} if record else None
    nsteps = int(round(sysd["event"]["simulated_duration_s"] / dt))
    T_old = 0.0  # motor torque gone; hydraulic torque evaluated below

    for step in range(1, nsteps + 1):
        Hn = [h.copy() for h in H]; Qn = [q.copy() for q in Q]
        for i in range(3):
            h, q = H[i], Q[i]
            CP = h[:-2] + B[i] * q[:-2] - R[i] * q[:-2] * np.abs(q[:-2])
            CM = h[2:] - B[i] * q[2:] + R[i] * q[2:] * np.abs(q[2:])
            Hn[i][1:-1] = 0.5 * (CP + CM)
            Qn[i][1:-1] = (CP - CM) / (2 * B[i])
        for i in range(2):  # junctions
            CP = H[i][-2] + B[i] * Q[i][-2] - R[i] * Q[i][-2] * abs(Q[i][-2])
            j = i + 1
            CM = H[j][1] - B[j] * Q[j][1] + R[j] * Q[j][1] * abs(Q[j][1])
            Hj = (CP / B[i] + CM / B[j]) / (1 / B[i] + 1 / B[j])
            Qj = (CP - Hj) / B[i]
            Hn[i][-1] = Hj; Qn[i][-1] = Qj; Hn[j][0] = Hj; Qn[j][0] = Qj
        CP = H[2][-2] + B[2] * Q[2][-2] - R[2] * Q[2][-2] * abs(Q[2][-2])
        Hn[2][-1] = Hres; Qn[2][-1] = (CP - Hres) / B[2]

        # upstream node: pump + check valve + vessel
        CM = H[0][1] - B[0] * Q[0][1] + R[0] * Q[0][1] * abs(Q[0][1])
        Tq_old = pump_T(Qp, a) if not cv_closed else TR * b0 * a * a

        def vessel_head(qv):
            vg = max(Vg + 0.5 * dt * (Qv + qv), 1e-6 * Vt)
            hg = Cgas / vg**n
            return hg - Hatm + zvb + (Vt - vg) / Av - (Cout if qv > 0 else Cin) * qv * abs(qv)

        def solve_open(a_new):
            qp, qv = Qp, Qv
            for _ in range(50):
                f1 = CM + B[0] * (qp + qv) - Hs - pump_h(qp, a_new)
                f2 = CM + B[0] * (qp + qv) - vessel_head(qv)
                e = 1e-7
                J11 = (CM + B[0] * (qp + e + qv) - Hs - pump_h(qp + e, a_new) - f1) / e
                J12 = B[0]
                J21 = B[0]
                J22 = (CM + B[0] * (qp + qv + e) - vessel_head(qv + e) - f2) / e
                det = J11 * J22 - J12 * J21
                dqp = (f1 * J22 - f2 * J12) / det
                dqv = (J11 * f2 - J21 * f1) / det
                lim = 0.05
                sc = max(1.0, abs(dqp) / lim, abs(dqv) / lim)
                dqp /= sc; dqv /= sc
                qp -= dqp; qv -= dqv
                if abs(dqp) + abs(dqv) < 1e-10:
                    break
            return qp, qv

        def solve_closed():
            qv = Qv
            for _ in range(50):
                f = CM + B[0] * qv - vessel_head(qv)
                e = 1e-7
                J = (CM + B[0] * (qv + e) - vessel_head(qv + e) - f) / e
                d = f / J
                d = max(-0.05, min(0.05, d))
                qv -= d
                if abs(d) < 1e-11:
                    break
            return qv

        if not cv_closed:
            a_new = a - dt * Tq_old / (I * wR)
            for _ in range(3):  # trapezoidal speed update
                qp, qv = solve_open(a_new)
                a_new = a - 0.5 * dt * (Tq_old + pump_T(max(qp, 0.0), a_new)) / (I * wR)
            # the check valve stays open only if the pump, at zero flow, can
            # still exceed the head it would see with the valve shut
            qv_shut = solve_closed()
            HP_shut = CM + B[0] * qv_shut
            if Hs + pump_h(0.0, a_new) <= HP_shut:
                cv_closed = True
            else:
                qp, qv = solve_open(a_new)
                if qp <= 0.0:
                    cv_closed = True
        if cv_closed:
            qp = 0.0
            qv = solve_closed()
            a_new = a - 0.5 * dt * (Tq_old + TR * b0 * a * a) / (I * wR)
        Vg = Vg + 0.5 * dt * (Qv + qv)
        Qv, Qp, a = qv, qp, max(a_new, 0.0)
        HP = CM + B[0] * (Qp + Qv)
        Hn[0][0] = HP; Qn[0][0] = Qp + Qv
        H, Q = Hn, Qn

        for i in range(3):
            pa = H[i] - z[i] + Hatm
            k = int(np.argmin(pa))
            if pa[k] < min_abs:
                min_abs = float(pa[k]); where_min = (step * dt, float(x0[i][k]))
            pg = H[i] - z[i]
            k = int(np.argmax(pg))
            if pg[k] > max_gauge:
                max_gauge = float(pg[k]); where_max = (step * dt, float(x0[i][k]))
        max_gas = max(max_gas, Vg / Vt)
        if not (math.isfinite(min_abs) and math.isfinite(max_gauge) and math.isfinite(Vg)) or Vg <= 0:
            min_abs, max_gauge, max_gas = -math.inf, math.inf, math.inf
            break
        if record and step % max(1, int(round(0.05 / dt))) == 0:
            for k, v in (("t", step * dt), ("Qp", Qp), ("Qv", Qv), ("HP", HP), ("alpha", a), ("Vg", Vg)):
                rec[k].append(v)
    out = {"min_absolute_pressure_head_m": min_abs, "max_gauge_pressure_head_m": max_gauge,
           "max_gas_fraction": float(max_gas), "steady_flow_m3_s": Q0, "wave_speed_adjustment": adj,
           "min_pressure_time_chainage": where_min, "max_pressure_time_chainage": where_max}
    if record:
        out["record"] = rec
    return out


def feasible(res, sysd, margin=0.0):
    lim = sysd["design_limits"]
    return (res["min_absolute_pressure_head_m"] >= lim["minimum_absolute_pressure_head_m"] + margin
            and res["max_gauge_pressure_head_m"] <= lim["maximum_gauge_pressure_head_m"] - margin
            and res["max_gas_fraction"] <= lim["maximum_gas_fraction_of_vessel"] - margin * 0.01)
