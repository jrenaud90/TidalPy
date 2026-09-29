"""The radial data-file reader and solid/liquid layer detection, on the bundled PREM and on hand-made profiles."""
import numpy as np
import pytest

from TidalPy.Structures.configs import data_file, worldpack


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
    assert radius[-1] == pytest.approx(6371.0e3)
    assert np.all(np.diff(radius) >= 0.0)


def test_prem_derived_moduli_match_formulas():
    arrays = data_file.load_radial_data(_prem_path())
    rho, vp, vs = arrays["density_kg_m3"], arrays["vp_m_s"], arrays["vs_m_s"]
    np.testing.assert_allclose(arrays["shear_modulus_pa"], rho * vs ** 2, rtol=1e-12)
    np.testing.assert_allclose(arrays["bulk_modulus_pa"], rho * (vp ** 2 - (4.0 / 3.0) * vs ** 2), rtol=1e-12)
    assert np.all(arrays["shear_modulus_pa"][vs == 0.0] == 0.0)


def test_prem_boundary_rows_keep_the_lower_layer_first():
    """At every duplicated radius the lower layer's row leads, though PREM lists the surface first."""
    arrays = data_file.load_radial_data(_prem_path())
    radius, density, shear = arrays["radius_m"], arrays["density_kg_m3"], arrays["shear_modulus_pa"]
    duplicated = np.flatnonzero(np.diff(radius) == 0.0)
    assert duplicated.size >= 2
    icb = duplicated[np.isclose(radius[duplicated], 1221.5e3)]
    assert icb.size == 1
    assert shear[icb[0]] > 0.0 and shear[icb[0] + 1] == 0.0
    assert density[icb[0]] > density[icb[0] + 1]
    cmb = duplicated[np.isclose(radius[duplicated], 3480.0e3)]
    assert cmb.size == 1
    assert shear[cmb[0]] == 0.0 and shear[cmb[0] + 1] > 0.0


def test_prem_center_first_copy_loads_identically(tmp_path):
    arrays = data_file.load_radial_data(_prem_path())
    # Radius, density, and the velocities: a headerless file reads its columns by position, and the fifth and sixth
    # positions are viscosities, not the quality factors PREM.csv adds.
    center_first = np.loadtxt(_prem_path(), delimiter=",", comments="#", usecols=range(4))[::-1]
    path = str(tmp_path / "center_first.csv")
    np.savetxt(path, center_first, delimiter=",")
    reloaded = data_file.load_radial_data(path)
    np.testing.assert_allclose(reloaded["radius_m"], arrays["radius_m"])
    np.testing.assert_allclose(reloaded["density_kg_m3"], arrays["density_kg_m3"])


def test_detect_three_layers_alternating():
    """PREM splits into inner core, outer core, and mantle plus crust, each ending on its own boundary row."""
    arrays = data_file.load_radial_data(_prem_path())
    radius, shear = arrays["radius_m"], arrays["shear_modulus_pa"]
    layers = data_file.detect_layer_boundaries(radius, shear)
    assert len(layers) == 3
    assert [is_solid for (_, _, is_solid) in layers] == [True, False, True]
    assert radius[layers[0][1]] == pytest.approx(1221.5e3) and shear[layers[0][1]] > 0.0
    assert radius[layers[1][1]] == pytest.approx(3480.0e3) and shear[layers[1][1]] == 0.0
    assert radius[layers[2][1]] == pytest.approx(6371.0e3)
    previous_outer = -1.0
    for start, end, _ in layers:
        assert radius[end] > radius[start]
        assert radius[end] > previous_outer
        previous_outer = radius[end]


@pytest.mark.parametrize("radius, shear, expected", [
    pytest.param([], [], [], id="empty"),
    pytest.param([0.0, 0.0, 1.0e6, 2.0e6], [0.0, 1.0e10, 1.0e10, 1.0e10], [(0, 3, True)], id="center-duplicate"),
    pytest.param(
        [0.0, 1.0e6, 1.0e6, 2.0e6],
        [1.0e10, 1.0e10, 0.0, 0.0],
        [(0, 1, True), (2, 3, False)],
        id="real-boundary"),
    # A stray liquid point at a duplicated radius joins the layer below it.
    pytest.param(
        [0.0, 1.0e6, 1.0e6, 2.0e6],
        [1.0e10, 0.0, 1.0e10, 1.0e10],
        [(0, 1, True), (2, 3, True)],
        id="stray-liquid-point"),
    pytest.param([0.0, 1.0e6, 2.0e6], [0.0, 0.0, 0.0], [(0, 2, False)], id="all-liquid"),
])
def test_no_layer_is_ever_zero_thickness(radius, shear, expected):
    radius = np.array(radius, dtype=np.float64)
    layers = data_file.detect_layer_boundaries(radius, np.array(shear, dtype=np.float64))
    assert layers == expected
    assert all(radius[end] > radius[start] for start, end, _ in layers)


# =====================================================================================================================
# The parser: names, units, layout
# =====================================================================================================================
def test_columns_are_found_by_name_in_any_order(tmp_path):
    path = _write(tmp_path, "named.tsv",
                  "Vs [m/s]\tdensity_kg_m3\tRadius_km\tVp [km/s]\n"
                  "0.0\t10000.0\t0.0\t10.0\n"
                  "3600.0\t5000.0\t3000.0\t11.0\n"
                  "3900.0\t4000.0\t6000.0\t12.0\n")
    arrays = data_file.load_radial_data(path)
    np.testing.assert_allclose(arrays["radius_m"], [0.0, 3.0e6, 6.0e6])
    np.testing.assert_allclose(arrays["vp_m_s"], [10.0e3, 11.0e3, 12.0e3])
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
    """The bundled file's prose line looks like a header (commas, a matching word) but is not read as one."""
    with open(_prem_path()) as handle:
        prose = [line for line in handle if line.startswith("#")][0]
    assert "Dziewonski" in prose
    assert data_file.load_radial_data(_prem_path())["radius_m"][-1] == pytest.approx(6371.0e3)


def test_headerless_file_is_read_positionally(tmp_path):
    path = _write(tmp_path, "bare.csv", "0,9000,10000,0\n1000,8000,10000,3000\n2000,7000,10000,3000\n")
    arrays = data_file.load_radial_data(path)
    np.testing.assert_allclose(arrays["radius_m"], [0.0, 1.0e6, 2.0e6])
    np.testing.assert_allclose(arrays["vs_m_s"], [0.0, 3000.0, 3000.0])


def test_radius_in_metres_is_left_alone(tmp_path):
    """An unlabelled radius is read as km below 100 km and as m above it."""
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
    # Vs above sqrt(3/4) Vp makes K = rho (Vp^2 - 4/3 Vs^2) negative.
    ({"radius_km": [0, 1], "density": [1000, 1000], "vp": [1e3, 1e3], "vs": [1e3, 1e3]},
     "non-positive bulk modulus"),
    ({"radius_km": [0, 1], "density": [1000, 1000], "vp": [1e4, 1e4], "vs": [0, 0], "eta": [1e19, -1.0]},
     "non-positive shear viscosity"),
    ({"radius_km": [0.0], "density": [1000.0], "vp": [1e4], "vs": [0.0]}, "at least 2"),
    ({"stuff": [0, 1], "things": [1, 1]}, "nothing this reader knows"),
    ({"radius_m": [1.0e6, 1.0e6], "density": [1.0e3, 1.0e3], "vp": [1.0e4, 1.0e4], "vs": [0.0, 0.0]},
     "spans no radius"),
    ({"radius_km": [0, 1], "r_km": [0, 1], "density": [1, 1], "vp": [1, 1], "vs": [0, 0]}, "more than once"),
])
def test_unusable_profiles_say_why(columns, message):
    with pytest.raises(ValueError, match=message):
        data_file.load_radial_data(columns)


@pytest.mark.parametrize("text, message", [
    pytest.param("1,2,3\n4,5,6\n", "at least 4 columns", id="too-few-columns"),
    pytest.param("1,2,3,4,5,6,7\n1,2,3,4,5,6,7\n", "at most 6", id="too-many-columns"),
    pytest.param("1,2,3,4\n4,5,6\n", "differing widths", id="ragged"),
    pytest.param("1,2,3,4\n4,5,six,7\n", "not numeric", id="non-numeric"),
])
def test_unusable_files_say_why(tmp_path, text, message):
    path = _write(tmp_path, "profile.csv", text)
    with pytest.raises(ValueError, match=message):
        data_file.load_radial_data(path)


def test_a_profile_must_be_a_path_or_a_mapping():
    with pytest.raises(TypeError, match="mapping"):
        data_file.load_radial_data(42)
