# Run from the bundle root.
import json, sys, time, numpy as np
sys.path.insert(0, "solution"); import solve as S
data = sys.argv[1] if len(sys.argv) > 1 else "."
truth = json.load(open(f"{data}/tests/expected.json"))
def score(out, tag):
    ef = [abs(out["blades"][k]["resonance_frequency_Hz"]/v["resonance_frequency_Hz"]-1) for k,v in truth["blades"].items()]
    ez = [out["blades"][k]["damping_ratio"]/v["damping_ratio"]-1 for k,v in truth["blades"].items()]
    ea = [out["blades"][k]["steady_resonant_amplitude_mm"]/v["steady_resonant_amplitude_mm"]-1 for k,v in truth["blades"].items()]
    ok = out["engine_order"]==truth["engine_order"] and max(ef)<=0.0025 and max(map(abs,ez))<=0.25 and abs(np.mean(ez))<=0.07 and max(map(abs,ea))<=0.15
    print(f"{tag:22s} EO={out['engine_order']} max|df/f|={max(ef):.5f} dz: max|.|={max(map(abs,ez)):.3f} mean={np.mean(ez):+.3f} | da: max|.|={max(map(abs,ea)):.3f} mean={np.mean(ea):+.3f} | PASS={ok}")
d = f"{data}/environment/data"
t=time.time(); score(S.analyse(d), "transient (reference)"); print("  time %.1fs" % (time.time()-t))
score(S.analyse(d, method="steady"), "steady-state SDOF")
score(S.analyse(d, use_southwell=False), "transient no-Southwell")
