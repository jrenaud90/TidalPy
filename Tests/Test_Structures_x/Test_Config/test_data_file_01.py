"""
Tests for the radial data-file reader and layer auto-detection
(``TidalPy.structures_x.configs.data_file``).

Covers the bundled ``PREM.csv`` (PREM with its 3 km ocean replaced by the upper crust): its MKS
arrays and derived moduli, and the solid/liquid layer detection that finds its three layers (inner
core, outer core, mantle plus crust). Covers the parser itself on the shapes real profiles come in:
columns found by name whatever their order, units taken from the column names or from the magnitude
of the radius, a header given as a comment or as a leading row or not at all, a depth column instead
of a radius, moduli instead of velocities, a mapping of arrays instead of a file, and the errors a
profile that cannot describe a planet produces.
"""

import numpy as np
import pytest

from TidalPy.structures_x.configs import data_file, worldpack


def _prem_path():
    return worldpack.resolve_data_file("PREM.csv")


def _write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text)
    return str(path)


# =====================================================================================================================
# The bundled PREM profile
# =====================================================================================================================
def test_load_prem_keys_and_shapes():
    arrays = data_file.load_radial_data(_prem_path())
    for key in ("radius_m", "density_kg_m3", "vp_m_s", "vs_m_s", "shear_modulus_pa", "bulk_modulus_pa"):
        assert arrays[key] is not None
        assert arrays[key].shape == arrays["radius_m"].shape
    assert arrays["radius_m"].size > 10
    # PREM names no viscosity, so the body it describes is elastic.
    assert arrays["shear_viscosity_pas"] is None
    assert arrays["bulk_viscosity_pas"] is None


def test_prem_radius_converted_to_metres_and_ascending():
    radius = data_file.load_radial_data(_prem_path())["radius_m"]
    assert radius[-1] == pytest.approx(6371.0e3)   # PREM's surface, given in km
    assert np.all(np.diff(radius) >= 0.0)


def test_prem_derived_moduli_match_formulas():
    arrays = data_file.load_radial_data(_prem_path())
    rho, vp, vs = arrays["density_kg_m3"], arrays["vp_m_s"], arrays["vs_m_s"]
    np.testing.assert_allclose(arrays["shear_modulus_pa"], rho * vs ** 2, rtol=1e-12)
    np.testing.assert_allclose(arrays["bulk_modulus_pa"], rho * (vp ** 2 - (4.0 / 3.0) * vs ** 2), rtol=1e-12)
    assert np.all(arrays["shear_modulus_pa"][vs == 0.0] == 0.0)


def test_prem_boundary_rows_keep_the_lower_layer_first():
    """PREM lists the surface first; after loading, the lower layer's row leads at every duplicated radius."""
    arrays = data_file.load_radial_data(_prem_path())
    radius, density, shear = arrays["radius_m"], arrays["density_kg_m3"], arrays["shear_modulus_pa"]
    duplicated = np.flatnonzero(np.diff(radius) == 0.0)
    assert duplicated.size >= 2
    # Inner-core boundary at 1221.5 km: the solid inner core (denser, non-zero shear) leads the liquid outer core.
    icb = duplicated[np.isclose(radius[duplicated], 1221.5e3)]
    assert icb.size == 1
    assert shear[icb[0]] > 0.0 and shear[icb[0] + 1] == 0.0
    assert density[icb[0]] > density[icb[0] + 1]
    # Core-mantle boundary at 3480 km: the liquid outer core leads the solid mantle.
    cmb = duplicated[np.isclose(radius[duplicated], 3480.0e3)]
    assert cmb.size == 1
    assert shear[cmb[0]] == 0.0 and shear[cmb[0] + 1] > 0.0


def test_prem_center_first_copy_loads_identically(tmp_path):
    arrays = data_file.load_radial_data(_prem_path())
    center_first = np.loadtxt(_prem_path(), delimiter=",", comments="#")[::-1]
    path = str(tmp_path / "center_first.csv")
    np.savetxt(path, center_first, delimiter=",")
    reloaded = data_file.load_radial_data(path)
    np.testing.assert_allclose(reloaded["radius_m"], arrays["radius_m"])
    np.testing.assert_allclose(reloaded["density_kg_m3"], arrays["density_kg_m3"])


def test_detect_three_layers_alternating():
    arrays = data_file.load_radial_data(_prem_path())
    radius, shear = arrays["radius_m"], arrays["shear_modulus_pa"]
    layers = data_file.detect_layer_boundaries(radius, shear)
    # Bundled PREM: inner core (solid), outer core (liquid), mantle plus crust (solid) to the surface.
    assert len(layers) == 3
    assert [is_solid for (_, _, is_solid) in layers] == [True, False, True]
    # Each layer ends on its own boundary row: the inner core on the dense solid row at 1221.5 km, the
    # outer core on the liquid row at 3480 km, the mantle at the surface.
    assert radius[layers[0][1]] == pytest.approx(1221.5e3) and shear[layers[0][1]] > 0.0
    assert radius[layers[1][1]] == pytest.approx(3480.0e3) and shear[layers[1][1]] == 0.0
    assert radius[layers[2][1]] == pytest.approx(6371.0e3)
    # Layers are non-degenerate (real radius span) and ordered inner to outer.
    previous_outer = -1.0
    for start, end, _ in layers:
        assert radius[end] > radius[start]
        assert radius[end] > previous_outer
        previous_outer = radius[end]


def test_detect_layer_boundaries_handles_empty():
    empty = np.empty(0, dtype=np.float64)
    assert data_file.detect_layer_boundaries(empty, empty) == []


@pytest.mark.parametrize("radius, shear, expected", [
    # Two rows at the very center with opposite solidity: one layer, not a zero-thickness one first.
    ([0.0, 0.0, 1.0e6, 2.0e6], [0.0, 1.0e10, 1.0e10, 1.0e10], [(0, 3, True)]),
    # A duplicated radius at a real boundary: each layer ends and starts on its own row.
    ([0.0, 1.0e6, 1.0e6, 2.0e6], [1.0e10, 1.0e10, 0.0, 0.0], [(0, 1, True), (2, 3, False)]),
    # A stray liquid point at a duplicated radius joins the layer below it.
    ([0.0, 1.0e6, 1.0e6, 2.0e6], [1.0e10, 0.0, 1.0e10, 1.0e10], [(0, 1, True), (2, 3, True)]),
    # A wholly liquid body is one liquid layer.
    ([0.0, 1.0e6, 2.0e6], [0.0, 0.0, 0.0], [(0, 2, False)]),
])
def test_no_layer_is_ever_zero_thickness(radius, shear, expected):
    radius = np.array(radius)
    layers = data_file.detect_layer_boundaries(radius, np.array(shear))
    assert layers == expected
    assert all(radius[end] > radius[start] for start, end, _ in layers)


def test_a_profile_must_span_a_radius_range():
    with pytest.raises(ValueError, match="spans no radius"):
        data_file.load_radial_data({"radius_m": [1.0e6, 1.0e6], "density": [1.0e3, 1.0e3],
                                    "vp": [1.0e4, 1.0e4], "vs": [0.0, 0.0]})


# =====================================================================================================================
# The parser: names, units, layout
# =====================================================================================================================
def test_columns_are_found_by_name_in_any_order(tmp_path):
    """A header row names the columns, so their order does not matter and they may state their units."""
    path = _write(tmp_path, "named.tsv",
                  "Vs [m/s]\tdensity_kg_m3\tRadius_km\tVp [km/s]\n"
                  "0.0\t10000.0\t0.0\t10.0\n"
                  "3600.0\t5000.0\t3000.0\t11.0\n"
                  "3900.0\t4000.0\t6000.0\t12.0\n")
    arrays = data_file.load_radial_data(path)
    np.testing.assert_allclose(arrays["radius_m"], [0.0, 3.0e6, 6.0e6])
    np.testing.assert_allclose(arrays["vp_m_s"], [10.0e3, 11.0e3, 12.0e3])   # km/s converted
    np.testing.assert_allclose(arrays["density_kg_m3"], [10000.0, 5000.0, 4000.0])


def test_header_may_be_the_last_comment_line(tmp_path):
    path = _write(tmp_path, "commented.csv",
                  "# a profile of something\n"
                  "# radius_m; rho; vp; vs; eta; eta_bulk\n"
                  "0.0; 9000.0; 10000.0; 0.0; 1e19; 1e20\n"
                  "1.0e6; 8000.0; 10000.0; 3000.0; 1e19; 1e20\n"
                  "2.0e6; 7000.0; 10000.0; 3000.0; 1e19; 1e20\n")
    arrays = data_file.load_radial_data(path)
    np.testing.assert_allclose(arrays["radius_m"], [0.0, 1.0e6, 2.0e6])
    np.testing.assert_allclose(arrays["shear_viscosity_pas"], 1.0e19)
    np.testing.assert_allclose(arrays["bulk_viscosity_pas"], 1.0e20)


def test_prose_comments_are_not_mistaken_for_a_header():
    """The bundled file's prose has the right number of commas and a matching word; it is still not a header."""
    with open(_prem_path()) as handle:
        prose = [line for line in handle if line.startswith("#")][0]
    assert "Dziewonski" in prose
    # It loads anyway, through the real header line below the prose.
    assert data_file.load_radial_data(_prem_path())["radius_m"][-1] == pytest.approx(6371.0e3)


def test_headerless_file_is_read_positionally(tmp_path):
    path = _write(tmp_path, "bare.csv", "0,9000,10000,0\n1000,8000,10000,3000\n2000,7000,10000,3000\n")
    arrays = data_file.load_radial_data(path)
    np.testing.assert_allclose(arrays["radius_m"], [0.0, 1.0e6, 2.0e6])   # km by magnitude
    np.testing.assert_allclose(arrays["vs_m_s"], [0.0, 3000.0, 3000.0])


def test_radius_in_metres_is_left_alone(tmp_path):
    """An unlabelled radius is read as km below 100 km and as m above it, ranges that cannot overlap."""
    path = _write(tmp_path, "metres.csv", "0,9000,10000,0\n1.0e6,8000,10000,3000\n2.0e6,7000,10000,3000\n")
    np.testing.assert_allclose(data_file.load_radial_data(path)["radius_m"], [0.0, 1.0e6, 2.0e6])


def test_depth_column_needs_the_world_radius(tmp_path):
    path = _write(tmp_path, "depth.csv",
                  "depth_km,density,shear_modulus,bulk_modulus\n"
                  "0,3000,6.0e10,1.0e11\n"
                  "1000,4000,7.0e10,2.0e11\n"
                  "2000,5000,8.0e10,3.0e11\n")
    arrays = data_file.load_radial_data(path, surface_radius=2.0e6)
    np.testing.assert_allclose(arrays["radius_m"], [0.0, 1.0e6, 2.0e6])
    # Given the moduli outright, the profile needs no velocities and reports none.
    np.testing.assert_allclose(arrays["shear_modulus_pa"], [8.0e10, 7.0e10, 6.0e10])
    assert arrays["vp_m_s"] is None and arrays["vs_m_s"] is None
    with pytest.raises(ValueError, match="radius_m"):
        data_file.load_radial_data(path)


def test_a_mapping_of_arrays_is_a_profile():
    """The Python route: no file, just the arrays."""
    arrays = data_file.load_radial_data({
        "radius_km": [0.0, 1000.0, 2000.0],
        "density":   [9000.0, 8000.0, 7000.0],
        "vp":        [10000.0, 10000.0, 10000.0],
        "vs":        [3000.0, 3000.0, 0.0],
    })
    np.testing.assert_allclose(arrays["radius_m"], [0.0, 1.0e6, 2.0e6])
    np.testing.assert_allclose(arrays["shear_modulus_pa"], [8.1e10, 7.2e10, 0.0])
    assert arrays["shear_viscosity_pas"] is None


# =====================================================================================================================
# The errors
# =====================================================================================================================
@pytest.mark.parametrize("columns, message", [
    ({"radius_km": [0, 1], "vp": [1, 1], "vs": [0, 0]}, "no density"),
    ({"radius_km": [0, 1], "density": [1, 1]}, "neither pair"),
    ({"density": [1, 1], "vp": [1, 1], "vs": [0, 0]}, "no radius"),
    ({"radius_km": [0, 1, 2], "density": [1, 1], "vp": [1, 1], "vs": [0, 0]}, "differing lengths"),
    ({"radius_km": [0, 1], "density": [1000, -1], "vp": [1, 1], "vs": [0, 0]}, "non-positive density"),
    # V_s above sqrt(3/4) V_p makes K = rho (Vp^2 - 4/3 Vs^2) negative: not a material.
    ({"radius_km": [0, 1], "density": [1000, 1000], "vp": [1e3, 1e3], "vs": [1e3, 1e3]},
     "non-positive bulk modulus"),
    ({"radius_km": [0, 1], "density": [1000, 1000], "vp": [1e4, 1e4], "vs": [0, 0], "eta": [1e19, -1.0]},
     "non-positive shear viscosity"),
    ({"radius_km": [0.0], "density": [1000.0], "vp": [1e4], "vs": [0.0]}, "at least 2"),
    ({"stuff": [0, 1], "things": [1, 1]}, "nothing this reader knows"),
])
def test_unusable_profiles_say_why(columns, message):
    with pytest.raises(ValueError, match=message):
        data_file.load_radial_data(columns)


def test_a_profile_must_be_a_path_or_a_mapping():
    with pytest.raises(TypeError, match="mapping"):
        data_file.load_radial_data(42)


def test_headerless_file_column_count_is_checked(tmp_path):
    too_few = _write(tmp_path, "few.csv", "1,2,3\n4,5,6\n")
    with pytest.raises(ValueError, match="at least 4 columns"):
        data_file.load_radial_data(too_few)
    too_many = _write(tmp_path, "many.csv", "1,2,3,4,5,6,7\n1,2,3,4,5,6,7\n")
    with pytest.raises(ValueError, match="at most 6"):
        data_file.load_radial_data(too_many)


def test_ragged_and_non_numeric_rows_are_rejected(tmp_path):
    ragged = _write(tmp_path, "ragged.csv", "1,2,3,4\n4,5,6\n")
    with pytest.raises(ValueError, match="differing widths"):
        data_file.load_radial_data(ragged)
    junk = _write(tmp_path, "junk.csv", "1,2,3,4\n4,5,six,7\n")
    with pytest.raises(ValueError, match="not numeric"):
        data_file.load_radial_data(junk)


def test_a_column_named_twice_is_rejected():
    with pytest.raises(ValueError, match="more than once"):
        data_file.load_radial_data({"radius_km": [0, 1], "r_km": [0, 1], "density": [1, 1],
                                    "vp": [1, 1], "vs": [0, 0]})
