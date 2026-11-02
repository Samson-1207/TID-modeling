"""Smoke tests: run with `pytest` from the repository root."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import config  # noqa: E402
from iri_msis_gw import iri_profiles, msis_profiles, gw_struct  # noqa: E402
from numerical_model import hooke_delta, I_deg  # noqa: E402

UT = np.array([0.0, 6.0, 12.0, 18.0])


@pytest.fixture(scope="module")
def model():
    N0, F2 = iri_profiles(UT, config.lat0, config.lon0, config.altF)
    ms = msis_profiles(UT, config.altF)
    gw = gw_struct(ms, config.altF)
    delta, uh = hooke_delta(N0[:, 0, :], gw, ms)
    return N0, F2, ms, gw, delta, uh


def test_background_shapes_and_values(model):
    N0, F2, ms, *_ = model
    assert N0.shape == (len(UT), 1, len(config.altF))
    assert np.all(N0 >= 1e6)
    assert np.all((F2["hm"] > 150) & (F2["hm"] < 450))
    assert np.all(ms["Tn"] > 100)


def test_gravity_wave_branch(model):
    gw = model[3]
    iref = np.argmin(np.abs(config.altF - config.z_ref))
    assert np.all(gw["m"].real < 0)          # downward phase, upward energy
    assert np.all(gw["m"].imag >= 0)         # damping with height
    assert np.allclose(gw["a"][:, iref], 1)  # normalised at z_ref
    assert np.all(gw["T_BV"] < config.period * 60)  # wave propagates vertically


def test_hooke_response_finite(model):
    delta, uh = model[4], model[5]
    assert delta.shape == (len(UT), len(config.altF))
    assert np.all(np.isfinite(delta))
    assert np.all(np.isfinite(uh))


def test_geomagnetic_field_tromso():
    assert 75 < I_deg < 80
