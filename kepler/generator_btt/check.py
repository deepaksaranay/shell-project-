import json, sys, time, numpy as np
sys.path.insert(0, "blade-tip-timing/solution"); import solve as S
truth = json.load(open("blade-tip-timing/tests/expected.json"))
def score(out, tag):
    ef = [abs(out["blades"][k]["resonance_frequency_Hz"]/v["resonance_frequency_Hz"]-1) for k,v in truth["blades"].items()]
    ez = [out["blades"][k]["damping_ratio"]/v["damping_ratio"]-1 for k,v in truth["blades"].items()]
    print(f"{tag:18s} EO={out['engine_order']} max|df/f|={max(ef):.5f} dz: max|.|={max(map(abs,ez)):.3f} mean={np.mean(ez):+.3f} rms={np.sqrt(np.mean(np.square(ez))):.3f}")
t=time.time(); score(S.analyse("blade-tip-timing/environment/data"), "reference"); print("time", time.time()-t)
score(S.analyse("blade-tip-timing/environment/data", use_southwell=False), "no_southwell")
