#!/usr/bin/env python3
"""Air-vessel sizing for the rising main in /app/data/system.json.

Independent transient model (written separately from any checking code):
  * MOC with a common time step and space-line interpolation (Courant <= 1)
    instead of wave-speed adjustment;
  * pump rotor speed integrated with Heun's method;
  * pump/check-valve/vessel node solved with scipy's fsolve;
  * check valve shuts for good once the pump's zero-flow head cannot exceed
    the head the line would have with the valve shut.
The volume is fixed just under the budget; the gas fraction and the inflow
orifice are scanned and the design with the most slack is written out.
"""
import json
import math
import os
import warnings

import numpy as np
from scipy.optimize import brentq, fsolve

SYSTEM_PATH = "/app/data/system.json"
DESIGN_PATH = "/app/output/surge_vessel_design.json"


class RisingMain:
    def __init__(self, cfg, dt):
        self.cfg, self.dt = cfg, dt
        f = cfg["fluid"]
        self.g, self.Ha = f["gravity_m_s2"], f["atmospheric_pressure_head_m"]
        p = cfg["pump"]
        self.QR, self.HR, self.TR = p["rated_flow_m3_s"], p["rated_head_m"], p["rated_torque_N_m"]
        self.Jw = p["polar_moment_of_inertia_kg_m2"] * p["rated_speed_rpm"] * math.pi / 30.0
        self.hc = [p["head_coefficients"][k] for k in ("h0", "h1", "h2")]
        self.tc = [p["torque_coefficients"][k] for k in ("b0", "b1", "b2")]
        self.Hsuc = cfg["suction_reservoir"]["water_level_m"]
        self.Hdis = cfg["discharge_reservoir"]["water_level_m"]
        prof = cfg["profile"]
        self.segs = []
        x_start = 0.0
        for pp in cfg["pipes"]:
            L, D, a, fr = pp["length_m"], pp["diameter_m"], pp["wave_speed_m_s"], pp["darcy_friction_factor"]
            n = max(1, int(math.floor(L / (a * dt))))  # Courant number a*dt/dx <= 1
            dx = L / n
            area = math.pi * D * D / 4
            x = x_start + dx * np.arange(n + 1)
            self.segs.append(dict(n=n, dx=dx, A=area, D=D, f=fr, a=a, theta=a * dt / dx,
                                  Bc=a / (self.g * area), Rc=fr * a * dt / (2 * self.g * D * area * area),
                                  z=np.interp(x, prof["chainage_m"], prof["elevation_m"])))
            x_start += L

    def head(self, q, al):
        v = q / self.QR
        return self.HR * (self.hc[0] * al * al + self.hc[1] * al * v + self.hc[2] * v * v)

    def torque(self, q, al):
        v = q / self.QR
        return self.TR * (self.tc[0] * al * al + self.tc[1] * al * v + self.tc[2] * v * v)

    def operating_point(self):
        g = self.g
        def gap(q):
            loss = sum(s["f"] * (s["dx"] * s["n"]) / s["D"] * (q / s["A"]) ** 2 / (2 * g) for s in self.segs)
            return self.Hsuc + self.head(q, 1.0) - self.Hdis - loss
        return brentq(gap, 1e-6, 3 * self.QR)

    def run(self, Vt, Vg0, d_in):
        cfg, dt, g, Ha = self.cfg, self.dt, self.g, self.Ha
        av = cfg["air_vessel"]
        nexp = av["gas_polytropic_exponent"]
        diam = (4 * Vt / (av["height_to_diameter_ratio"] * math.pi)) ** (1 / 3)
        Asec = math.pi * diam * diam / 4
        zb = av["bottom_elevation_m"]
        Aout = math.pi * av["outflow_connection"]["diameter_m"] ** 2 / 4
        k_out = av["outflow_connection"]["loss_coefficient"] / (2 * g * Aout * Aout)
        k_in = av["inflow_orifice"]["loss_coefficient"] / (2 * g * (math.pi * d_in * d_in / 4) ** 2)

        q0 = self.operating_point()
        Hs, Qs = [], []
        h_up = self.Hsuc + self.head(q0, 1.0)
        for s in self.segs:
            grad = s["f"] / s["D"] * (q0 / s["A"]) ** 2 / (2 * g)
            Hs.append(h_up - grad * s["dx"] * np.arange(s["n"] + 1))
            Qs.append(np.full(s["n"] + 1, q0))
            h_up = Hs[-1][-1]

        level0 = zb + (Vt - Vg0) / Asec
        Kgas = (Hs[0][0] + Ha - level0) * Vg0 ** nexp
        vg, qv, qp, al, shut = Vg0, 0.0, q0, 1.0, False
        lo = min(float(np.min(H - s["z"])) for H, s in zip(Hs, self.segs)) + Ha
        hi = max(float(np.max(H - s["z"])) for H, s in zip(Hs, self.segs))
        gmax = vg / Vt

        def interp_side(H, Q, s, toward):
            """Values at the foot of a characteristic inside the reach."""
            th = s["theta"]
            if toward == "left":   # C+ foot between j-1 and j
                Hf = H[1:] - th * (H[1:] - H[:-1]); Qf = Q[1:] - th * (Q[1:] - Q[:-1])
            else:                  # C- foot between j and j+1
                Hf = H[:-1] - th * (H[:-1] - H[1:]); Qf = Q[:-1] - th * (Q[:-1] - Q[1:])
            return Hf, Qf

        steps = int(round(cfg["event"]["simulated_duration_s"] / dt))
        for _ in range(steps):
            Cp, Cm = [], []
            for H, Q, s in zip(Hs, Qs, self.segs):
                Hr, Qr = interp_side(H, Q, s, "left")
                Hsd, Qsd = interp_side(H, Q, s, "right")
                Cp.append(Hr + s["Bc"] * Qr - s["Rc"] * Qr * np.abs(Qr))     # valid at nodes 1..n
                Cm.append(Hsd - s["Bc"] * Qsd + s["Rc"] * Qsd * np.abs(Qsd))  # valid at nodes 0..n-1
            newH, newQ = [], []
            for k, s in enumerate(self.segs):
                H = np.empty(s["n"] + 1); Q = np.empty(s["n"] + 1)
                H[1:-1] = 0.5 * (Cp[k][:-1] + Cm[k][1:])
                Q[1:-1] = (Cp[k][:-1] - Cm[k][1:]) / (2 * s["Bc"])
                newH.append(H); newQ.append(Q)
            for k in range(len(self.segs) - 1):
                b1, b2 = self.segs[k]["Bc"], self.segs[k + 1]["Bc"]
                cp, cm = Cp[k][-1], Cm[k + 1][0]
                Hj = (cp / b1 + cm / b2) / (1 / b1 + 1 / b2)
                newH[k][-1] = newH[k + 1][0] = Hj
                newQ[k][-1] = newQ[k + 1][0] = (cp - Hj) / b1
            sl = self.segs[-1]
            newH[-1][-1] = self.Hdis
            newQ[-1][-1] = (Cp[-1][-1] - self.Hdis) / sl["Bc"]

            cm0, b0 = Cm[0][0], self.segs[0]["Bc"]

            def tee_head(qv_new):
                gas = max(vg + 0.5 * dt * (qv + qv_new), 1e-9)
                kk = k_out if qv_new > 0 else k_in
                return Kgas / gas ** nexp - Ha + zb + (Vt - gas) / Asec - kk * qv_new * abs(qv_new)

            def shut_flow():
                return fsolve(lambda u: cm0 + b0 * u[0] - tee_head(u[0]), [qv], xtol=1e-12)[0]

            if not shut:
                tq0 = self.torque(qp, al)
                al_pred = al - dt * tq0 / self.Jw
                qv_s = shut_flow()
                if self.Hsuc + self.head(0.0, al_pred) <= cm0 + b0 * qv_s:
                    shut = True
                else:
                    def eqs(u, al_n):
                        a_, b_ = u
                        return [cm0 + b0 * (a_ + b_) - self.Hsuc - self.head(a_, al_n),
                                cm0 + b0 * (a_ + b_) - tee_head(b_)]
                    qp_n, qv_n = fsolve(eqs, [qp, qv], args=(al_pred,), xtol=1e-12)
                    al_new = al - 0.5 * dt * (tq0 + self.torque(max(qp_n, 0.0), al_pred)) / self.Jw
                    qp_n, qv_n = fsolve(eqs, [qp_n, qv_n], args=(al_new,), xtol=1e-12)
                    if qp_n <= 0:
                        shut = True
                    else:
                        vg += 0.5 * dt * (qv + qv_n)
                        qp, qv, al = qp_n, qv_n, al_new
            if shut:
                qv_n = shut_flow()
                vg += 0.5 * dt * (qv + qv_n)
                al = max(al - dt * self.TR * self.tc[0] * al * al / self.Jw, 0.0)
                qp, qv = 0.0, qv_n
            newH[0][0] = cm0 + b0 * (qp + qv)
            newQ[0][0] = qp + qv
            Hs, Qs = newH, newQ
            for H, s in zip(Hs, self.segs):
                lo = min(lo, float(np.min(H - s["z"])) + Ha)
                hi = max(hi, float(np.max(H - s["z"])))
            gmax = max(gmax, vg / Vt)
            if not (math.isfinite(lo) and math.isfinite(hi) and vg > 0):
                return -math.inf, math.inf, math.inf
        return lo, hi, gmax


def main():
    with open(SYSTEM_PATH) as fh:
        cfg = json.load(fh)
    lim = cfg["design_limits"]
    model = RisingMain(cfg, dt=0.01)
    Vt = round(cfg["design_budget"]["maximum_vessel_total_volume_m3"] - 0.1, 2)
    best = None
    for d_in in (0.08, 0.12):
        for frac in np.arange(0.40, 0.49, 0.005):
            lo, hi, gmax = model.run(Vt, Vt * frac, d_in)
            slack = min(lo - lim["minimum_absolute_pressure_head_m"],
                        lim["maximum_gauge_pressure_head_m"] - hi,
                        25.0 * (lim["maximum_gas_fraction_of_vessel"] - gmax))
            if best is None or slack > best[0]:
                best = (slack, frac, d_in, lo, hi, gmax)
    slack, frac, d_in, lo, hi, gmax = best
    design = {"vessel_total_volume_m3": Vt,
              "initial_gas_volume_m3": round(Vt * frac, 3),
              "inflow_orifice_diameter_m": d_in}
    os.makedirs(os.path.dirname(DESIGN_PATH), exist_ok=True)
    with open(DESIGN_PATH, "w") as fh:
        json.dump(design, fh, indent=2)
    print(json.dumps(design), f"slack {slack:.3f} | min abs head {lo:.2f} m, max gauge head {hi:.2f} m, gas {gmax:.3f}")


if __name__ == "__main__":
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    main()
