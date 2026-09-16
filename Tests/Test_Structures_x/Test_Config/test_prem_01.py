"""
Tests for the Cython PREM-like data loader and layer auto-detection
(``TidalPy.structures_x.configs.prem``).

Loads the bundled ``PREM.csv`` (PREM with its 3 km ocean replaced by the upper crust) and checks
that the MKS arrays and derived moduli are correct, that a surface-first file keeps the lower
layer's row first at a duplicated boundary radius, and that the solid/liquid layer
auto-detection finds the three layers (inner core, outer core, mantle plus crust) from the shear
modulus.
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
    radius = arrays["radius_m"]
    # PREM surface is ~6.371e6 m.
    assert radius[-1] == pytest.approx(6371.0e3)
    assert np.all(np.diff(radius) >= 0.0)


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


def test_boundary_rows_keep_the_lower_layer_first():
    """PREM lists the surface first; after loading, the lower layer's row leads at every duplicated radius."""
    arrays = prem.load_prem_arrays(_prem_path())
    radius = arrays["radius_m"]
    density = arrays["density_kg_m3"]
    shear = arrays["shear_modulus_pa"]
    duplicated = np.flatnonzero(np.diff(radius) == 0.0)
    assert duplicated.size >= 2
    # Inner-core boundary at 1221.5 km: the solid inner core (denser, nonzero shear) precedes the liquid outer core.
    icb = duplicated[np.isclose(radius[duplicated], 1221.5e3)]
    assert icb.size == 1
    assert shear[icb[0]] > 0.0 and shear[icb[0] + 1] == 0.0
    assert density[icb[0]] > density[icb[0] + 1]
    # Core-mantle boundary at 3480 km: the liquid outer core precedes the solid mantle.
    cmb = duplicated[np.isclose(radius[duplicated], 3480.0e3)]
    assert cmb.size == 1
    assert shear[cmb[0]] == 0.0 and shear[cmb[0] + 1] > 0.0

    # A center-first copy of the file loads to the same arrays.
    center_first = np.loadtxt(_prem_path(), delimiter=",", comments="#")[::-1]
    path = _prem_path() + ".center_first.tmp"
    try:
        np.savetxt(path, center_first, delimiter=",")
        reloaded = prem.load_prem_arrays(path)
    finally:
        import os
        os.remove(path)
    np.testing.assert_allclose(reloaded["radius_m"], radius)
    np.testing.assert_allclose(reloaded["density_kg_m3"], density)


def test_detect_three_layers_alternating():
    arrays = prem.load_prem_arrays(_prem_path())
    radius = np.ascontiguousarray(arrays["radius_m"], dtype=np.float64)
    shear = np.ascontiguousarray(arrays["shear_modulus_pa"], dtype=np.float64)
    layers = prem.detect_layer_boundaries(radius, shear)
    # Bundled PREM: inner core (solid), outer core (liquid), mantle plus crust (solid) to the surface.
    assert len(layers) == 3
    solidity = [is_solid for (_, _, is_solid) in layers]
    assert solidity == [True, False, True]
    # Each layer ends on its own boundary row: the inner core on the dense solid row at 1221.5 km, the outer
    # core on the liquid row at 3480 km, the mantle at the surface.
    assert radius[layers[0][1]] == pytest.approx(1221.5e3) and shear[layers[0][1]] > 0.0
    assert radius[layers[1][1]] == pytest.approx(3480.0e3) and shear[layers[1][1]] == 0.0
    assert radius[layers[2][1]] == pytest.approx(6371.0e3)
    # Layers are non-degenerate (real radius span) and ordered inner-to-outer.
    prev_outer = -1.0
    for start, end, _ in layers:
        assert radius[end] > radius[start]
        assert radius[end] > prev_outer
        prev_outer = radius[end]


def test_no_data_file_handles_empty():
    # Detection on empty arrays returns no layers (defensive).
    empty = np.empty(0, dtype=np.float64)
    assert prem.detect_layer_boundaries(empty, empty) == []
