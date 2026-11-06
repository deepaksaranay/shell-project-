# Alternative correct method: time-domain least-squares sinusoid fit (no window),
# with 1X..3X runout regressors at the keyphasor speed; fitted (not raw) support model.
import json, sys, csv, os, numpy as np
sys.path.insert(0, "bearing-impedance-identification/solution")
import solve as S
D = "bearing-impedance-identification/environment/data"
truth = json.load(open("bearing-impedance-identification/tests/expected.json"))
w = 2*np.pi*S.TONES
def amps(r):
    t = r["t"]; fr = S.shaft_speed_hz(r)
    freqs = list(S.TONES) + ([h*fr for h in (1,2,3)] if fr else [])
    B = np.column_stack([np.ones_like(t)] + [f(2*np.pi*fq*t) for fq in freqs for f in (np.cos, np.sin)])
    dx, dy = S.probes_to_xy(r, {})
    out = {}
    for k, s in [("fx",r["fx"]),("fy",r["fy"]),("ax",r["ax"]),("ay",r["ay"]),("dx",dx),("dy",dy)]:
        c, *_ = np.linalg.lstsq(B, s, rcond=None)
        n = len(S.TONES); out[k] = c[1:2*n+1:2] - 1j*c[2:2*n+2:2]
    return out
bx, by = amps(S.load(f"{D}/baseline_dry_x.csv")), amps(S.load(f"{D}/baseline_dry_y.csv"))
F0, A0 = S.mat(bx,by,"fx","fy"), S.mat(bx,by,"ax","ay"); X0 = -A0/(w**2)[:,None,None]
Hraw = np.swapaxes(np.linalg.solve(np.swapaxes(X0,1,2), np.swapaxes(F0-S.MS*A0,1,2)),1,2)
Ks, Cs, Ms = S.fit_kcm(Hraw, w, np.ones(len(w),bool))
Hs = Ks[None] - (w**2)[:,None,None]*Ms[None] + 1j*w[:,None,None]*Cs[None]
worst = [0,0]
for row in csv.DictReader(open(f"{D}/test_log.csv")):
    sx, sy = amps(S.load(f"{D}/{row['x_shaker_file']}")), amps(S.load(f"{D}/{row['y_shaker_file']}"))
    F, A, Dd = S.mat(sx,sy,"fx","fy"), S.mat(sx,sy,"ax","ay"), S.mat(sx,sy,"dx","dy")
    rhs = F - S.MS*A - Hs @ (-A/(w**2)[:,None,None])
    Hb = np.swapaxes(np.linalg.solve(np.swapaxes(Dd,1,2), np.swapaxes(rhs,1,2)),1,2)
    K, C, M = S.fit_kcm(Hb, w, np.ones(len(w),bool))
    v = truth[row["condition"]]; Kt, Ct = np.array(v["K"]), np.array(v["C"])
    worst[0] = max(worst[0], np.abs(K-Kt).max()/max(abs(Kt[0,0]),abs(Kt[1,1])))
    worst[1] = max(worst[1], np.abs(C-Ct).max()/max(abs(Ct[0,0]),abs(Ct[1,1])))
print("alt LS-sine method worst K, C:", worst)
