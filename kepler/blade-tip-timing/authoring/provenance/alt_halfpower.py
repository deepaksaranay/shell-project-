# Run from the bundle root.
# Alternative correct method: half-power bandwidth on smoothed |Q - c|, with the
# Southwell detuning-rate correction  zeta = (EO - S/EO) dOmega / (2 EO Omega_res).
import sys, json, numpy as np
sys.path.insert(0, "solution"); import solve as S
T = json.load(open("tests/expected.json"))
opr, toas = S.load("environment/data")
D, Om = S.deflections(opr, toas); D = S.zero_static(D, Om)
eo, Q = S.choose_eo(D); Q = S.clean_Q(Q)
ez, ef = [], []
for b in range(S.B):
    q = Q[:, b]; m = np.isfinite(q); o, q = Om[m], q[m]
    sm = np.convolve(np.abs(q), np.ones(15)/15, "same"); pk = o[np.argmax(sm)]
    far = np.abs(o - pk) > 4.0
    c = np.median(q[far].real) + 1j*np.median(q[far].imag)
    a = np.abs(np.convolve(q - c, np.ones(31)/31, "same"))
    k = np.argmax(a); half = a[k]/np.sqrt(2)
    lo = k
    while a[lo] > half: lo -= 1
    hi = k
    while a[hi] > half: hi += 1
    o1 = np.interp(half, [a[lo], a[lo+1]], [o[lo], o[lo+1]])
    o2 = np.interp(half, [a[hi], a[hi-1]], [o[hi], o[hi-1]])
    om_res = 0.5*(o1+o2)
    z = (eo - S.S/eo)*(o2-o1)/(2*eo*om_res)
    t = T["blades"][str(b+1)]
    ez.append(z/t["damping_ratio"]-1); ef.append(abs(eo*om_res/t["resonance_frequency_Hz"]-1))
print(f"half-power EO={eo} max|df/f|={max(ef):.5f} dz max={max(map(abs,ez)):.3f} mean={np.mean(ez):+.3f} rms={np.sqrt(np.mean(np.square(ez))):.3f}")
