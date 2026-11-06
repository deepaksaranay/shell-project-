import json
import math
import os
import pathlib

import pytest

OUT = pathlib.Path(os.environ.get("BEARING_OUTPUT", "/app/output/bearing_coefficients.json"))
TRUTH = json.loads((pathlib.Path(__file__).parent / "expected.json").read_text())

DIRECT_REL_TOL = 0.03   # direct terms: 3 % of the true value
CROSS_REL_TOL = 0.03    # cross terms: 3 % of the mean direct magnitude
# Added mass is the weakest-identified term (it enters only via w^2), so it gets a looser band.
MASS_DIRECT_REL_TOL = 0.10
MASS_CROSS_ABS_TOL = 0.5  # kg


@pytest.fixture(scope="module")
def result():
    assert OUT.is_file(), f"missing output file {OUT}"
    data = json.loads(OUT.read_text())
    assert isinstance(data, dict)
    return data


def test_keys_present_and_numeric(result):
    for p in "KCM":
        for ij in ("xx", "xy", "yx", "yy"):
            key = p + ij
            assert key in result, f"missing {key}"
            v = result[key]
            assert isinstance(v, (int, float)) and not isinstance(v, bool), key
            assert math.isfinite(v), key


def test_excitation_frequencies(result):
    freqs = result.get("frequencies_hz")
    assert isinstance(freqs, list) and freqs, "frequencies_hz must be a non-empty list"
    got = sorted(round(float(f), 3) for f in freqs)
    want = sorted(round(float(f), 3) for f in TRUTH["frequencies_hz"])
    assert got == want, f"excitation frequencies {got} != {want}"


@pytest.mark.parametrize("p", ["K", "C"])
def test_direct_coefficients(result, p):
    for ij in ("xx", "yy"):
        key = p + ij
        t = TRUTH[key]
        assert abs(result[key] - t) <= DIRECT_REL_TOL * abs(t), f"{key}: got {result[key]:.6g}, expected {t:.6g}"


@pytest.mark.parametrize("p", ["K", "C"])
def test_cross_coupled_coefficients(result, p):
    scale = 0.5 * (abs(TRUTH[p + "xx"]) + abs(TRUTH[p + "yy"]))
    for ij in ("xy", "yx"):
        key = p + ij
        t = TRUTH[key]
        assert abs(result[key] - t) <= CROSS_REL_TOL * scale, f"{key}: got {result[key]:.6g}, expected {t:.6g}"


def test_added_mass(result):
    for ij in ("xx", "yy"):
        key = "M" + ij
        t = TRUTH[key]
        assert abs(result[key] - t) <= MASS_DIRECT_REL_TOL * abs(t), f"{key}: got {result[key]:.6g}, expected {t:.6g}"
    for ij in ("xy", "yx"):
        key = "M" + ij
        t = TRUTH[key]
        assert abs(result[key] - t) <= MASS_CROSS_ABS_TOL, f"{key}: got {result[key]:.6g}, expected {t:.6g}"
