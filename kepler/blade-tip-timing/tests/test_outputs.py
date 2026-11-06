import json
import math
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import pytest

PROGRAM = Path("/app/btt_analyze.py")
HERE = Path(__file__).parent
# heldout/ holds a copy of the shipped example run plus the unseen runs.
DATASETS = {d.name: d for d in sorted((HERE / "heldout").iterdir()) if d.is_dir()}
TRUTH = {n: json.loads((HERE / "truth" / f"{n}.json").read_text()) for n in DATASETS}

RUN_TIMEOUT_S = 600
F_TOL = 0.0025
Z_TOL = 0.25
Z_MEAN_TOL = 0.07
A_TOL = 0.15
RUN_AS = int(os.environ.get("BTT_RUN_UID", "65534"))  # nobody


def _run(name):
    """Run the program on one dataset as an unprivileged user; return parsed output."""
    work = Path(tempfile.mkdtemp(prefix=f"btt_{name}_"))
    data = work / "data"
    shutil.copytree(DATASETS[name], data)
    out = work / "out.json"
    for p in [work, data, *data.iterdir()]:
        os.chmod(p, 0o777 if p.is_dir() else 0o644)
    kw = {}
    if os.geteuid() == 0:
        kw = dict(user=RUN_AS, group=RUN_AS, extra_groups=[])
    proc = subprocess.run(
        ["python3", str(PROGRAM), str(data), str(out)],
        cwd=work, capture_output=True, text=True, timeout=RUN_TIMEOUT_S,
        env={"PATH": os.environ.get("PATH", "/usr/local/bin:/usr/bin:/bin"), "HOME": str(work),
             "OMP_NUM_THREADS": "2", "OPENBLAS_NUM_THREADS": "2", "MKL_NUM_THREADS": "2"},
        **kw,
    )
    assert proc.returncode == 0, f"{name}: exit {proc.returncode}\n{proc.stderr[-2000:]}"
    assert out.exists(), f"{name}: no output written"
    res = json.loads(out.read_text())
    assert isinstance(res, dict), f"{name}: top level must be an object"
    return res


@pytest.fixture(scope="session")
def results():
    assert PROGRAM.exists(), f"{PROGRAM} not found"
    cache = {}
    for name in DATASETS:
        try:
            cache[name] = _run(name)
        except Exception as e:  # noqa: BLE001 - recorded and reported per dataset
            cache[name] = e
    return cache


def _get(results, name):
    r = results[name]
    if isinstance(r, Exception):
        pytest.fail(f"{name}: program failed: {r}")
    return r


def _values(res, truth, key, name):
    blades = res.get("blades")
    assert isinstance(blades, dict), f"{name}: blades must be an object"
    vals = {}
    for b in truth["blades"]:
        assert b in blades and isinstance(blades[b], dict), f"{name}: blade {b} missing"
        v = blades[b].get(key)
        assert isinstance(v, (int, float)) and not isinstance(v, bool), f"{name}: blade {b} {key} not a number"
        v = float(v)
        assert math.isfinite(v) and v > 0, f"{name}: blade {b} {key} not finite/positive"
        vals[b] = v
    return vals


@pytest.mark.parametrize("name", list(DATASETS))
def test_engine_order(results, name):
    res = _get(results, name)
    eo = res.get("engine_order")
    assert isinstance(eo, int) and not isinstance(eo, bool), f"{name}: engine_order must be an integer"
    assert eo == TRUTH[name]["engine_order"], f"{name}: engine order {eo}"


@pytest.mark.parametrize("name", list(DATASETS))
def test_resonance_frequency(results, name):
    t = TRUTH[name]
    v = _values(_get(results, name), t, "resonance_frequency_Hz", name)
    bad = {b: round(v[b] / t["blades"][b]["resonance_frequency_Hz"] - 1, 5) for b in v
           if abs(v[b] / t["blades"][b]["resonance_frequency_Hz"] - 1) > F_TOL}
    assert not bad, f"{name}: frequency out of limit for blades {bad}"


@pytest.mark.parametrize("name", list(DATASETS))
def test_damping_ratio(results, name):
    t = TRUTH[name]
    v = _values(_get(results, name), t, "damping_ratio", name)
    rel = {b: v[b] / t["blades"][b]["damping_ratio"] - 1 for b in v}
    bad = {b: round(r, 3) for b, r in rel.items() if abs(r) > Z_TOL}
    assert not bad, f"{name}: damping out of limit for blades {bad}"
    mean = sum(rel.values()) / len(rel)
    assert abs(mean) <= Z_MEAN_TOL, f"{name}: mean relative damping error {mean:+.3f}"


@pytest.mark.parametrize("name", list(DATASETS))
def test_steady_resonant_amplitude(results, name):
    t = TRUTH[name]
    v = _values(_get(results, name), t, "steady_resonant_amplitude_mm", name)
    bad = {b: round(v[b] / t["blades"][b]["steady_resonant_amplitude_mm"] - 1, 3) for b in v
           if abs(v[b] / t["blades"][b]["steady_resonant_amplitude_mm"] - 1) > A_TOL}
    assert not bad, f"{name}: amplitude out of limit for blades {bad}"
