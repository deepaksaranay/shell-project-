# Run from the bundle root:  python3 authoring/provenance/check.py [bundle_root]
# Scores the reference, a steady-state SDOF fit (measured peak amplitude) and a
# no-Southwell variant on every dataset against the stated limits.
import json, os, sys, time, warnings
import numpy as np
warnings.filterwarnings("ignore")
sys.path.insert(0, "solution"); import solve as S

root = sys.argv[1] if len(sys.argv) > 1 else "."
DATASETS = {}
for n in sorted(os.listdir(f"{root}/tests/heldout")):
    DATASETS[n] = f"{root}/tests/heldout/{n}"
LIM = dict(f=0.0025, z=0.25, zmean=0.07, a=0.15)


def score(out, truth):
    bl = list(truth["blades"])
    g = lambda b, k: out["blades"][b][k] / truth["blades"][b][k] - 1
    ef = max(abs(g(b, "resonance_frequency_Hz")) for b in bl)
    ez = [g(b, "damping_ratio") for b in bl]
    ea = [g(b, "steady_resonant_amplitude_mm") for b in bl]
    ok = (out["engine_order"] == truth["engine_order"] and ef <= LIM["f"] and max(map(abs, ez)) <= LIM["z"]
          and abs(np.mean(ez)) <= LIM["zmean"] and max(map(abs, ea)) <= LIM["a"])
    return ok, (f"EO={out['engine_order']}/{truth['engine_order']} df={ef:.5f} dz max={max(map(abs, ez)):.3f} "
                f"mean={np.mean(ez):+.3f} da max={max(map(abs, ea)):.3f} mean={np.mean(ea):+.3f}")


for name, d in DATASETS.items():
    truth = json.load(open(f"{root}/tests/truth/{name}.json"))
    for tag, kw in [("reference", {}), ("steady+peak", {"method": "steady"}), ("no-Southwell", {"use_southwell": False})]:
        t = time.time()
        try:
            ok, msg = score(S.analyse(d, **kw), truth)
        except Exception as e:  # noqa: BLE001
            ok, msg = False, f"crashed: {type(e).__name__}"
        print(f"{name:8s} {tag:13s} {'PASS' if ok else 'fail'}  {msg}  ({time.time()-t:.1f}s)", flush=True)
