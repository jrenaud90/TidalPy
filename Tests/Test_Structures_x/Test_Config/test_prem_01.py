"""
Tests for the Cython PREM-like data loader and layer auto-detection
(``TidalPy.structures_x.configs.prem``).

Loads the bundled ``PREM.csv`` and checks that the MKS arrays and derived moduli are
correct and that the solid/liquid layer auto-detection finds Earth's four layers
(inner core, outer core, mantle, ocean) from the shear modulus.
"""

import numpy as np
import pytest

from TidalPy.structures_x.configs import prem, worldpack


def _prem_path():
    return worldpack.resolve_data_file("PREM.csv")


def test_load_prem_arrays_keys_and_shapes():
    arrays = prem.load_prem_arrays(_prem_path())
    for key in ("radius_m", "density_kg_m3", "vp_m_s", "vs_m_s",
                "shear_modulus_pa", "bulk_modulus_pa"):
        assert key in arrays
        assert arrays[key] is not None
    n = arrays["radius_m"].shape[0]
    assert n > 10
    for key in ("density_kg_m3", "vp_m_s", "vs_m_s", "shear_modulus_pa", "bulk_modulus_pa"):
        assert arrays[key].shape[0] == n


def test_radius_converted_to_metres_and_ascending():
    arrays = prem.load_prem_arrays(_prem_path())
    radius_m = arrays["radius_m"]
    # PREM surface is ~6.371e6 m.
    assert radius_m[-1] == pytest.approx(6371.0e3)
    assert np.all(np.diff(radius_m) >= 0.0)


def test_derived_moduli_match_formulas():
    arrays = prem.load_prem_arrays(_prem_path())
    rho    = arrays["density_kg_m3"]
    vp     = arrays["vp_m_s"]
    vs     = arrays["vs_m_s"]
    expected_shear = rho * vs ** 2
    expected_bulk  = rho * (vp ** 2 - (4.0 / 3.0) * vs ** 2)
    np.testing.assert_allclose(arrays["shear_modulus_pa"], expected_shear, rtol=1e-12)
    np.testing.assert_allclose(arrays["bulk_modulus_pa"], expected_bulk, rtol=1e-12)
    # Liquid regions (Vs == 0) have zero shear modulus.
    assert np.all(arrays["shear_modulus_pa"][vs == 0.0] == 0.0)


def test_detect_four_layers_alternating():
    arrays = prem.load_prem_arrays(_prem_path())
    radius_m = np.ascontiguousarray(arrays["radius_m"], dtype=np.float64)
    shear = np.ascontiguousarray(arrays["shear_modulus_pa"], dtype=np.float64)
    layers = prem.detect_layer_boundaries(radius_m, shear)
    # PREM Earth: inner core (solid), outer core (liquid), mantle (solid), ocean (liquid).
    assert len(layers) == 4
    solidity = [is_solid for (_, _, is_solid) in layers]
    assert solidity == [True, False, True, False]
    # Layers are non-degenerate (real radius span) and ordered inner-to-outer.
    prev_outer = -1.0
    for start, end, _ in layers:
        assert radius_m[end] > radius_m[start]
        assert radius_m[end] > prev_outer
        prev_outer = radius_m[end]


def test_no_data_file_handles_empty():
    # Detection on empty arrays returns no layers (defensive).
    empty = np.empty(0, dtype=np.float64)
    assert prem.detect_layer_boundaries(empty, empty) == []
