"""Tests for the native `_x` radial-solver input builders (TidalPy.RadialSolver_x.build_inputs).

The builders mirror the classic ``TidalPy.RadialSolver.helpers`` builders but evaluate the complex
moduli with ``rheology_x`` models, accept a single model for every layer, and validate their inputs.
"""
import numpy as np
import pytest

from TidalPy.initialize import build_logging_x_config
from TidalPy.RadialSolver import radial_solver as radial_solver_old
from TidalPy.RadialSolver.helpers import build_rs_input_from_data as classic_from_data
from TidalPy.RadialSolver.helpers import build_rs_input_homogeneous_layers as classic_homogeneous
from TidalPy.RadialSolver_x import (
    PlanetBuildData,
    build_rs_input_from_data,
    build_rs_input_homogeneous_layers,
    radial_solver,
)
from TidalPy.rheology.models import Elastic as ClassicElastic
from TidalPy.rheology.models import Maxwell as ClassicMaxwell
from TidalPy.rheology_x import Andrade, Elastic, Maxwell
from TidalPy.Utilities_x.logging_x.logger import flush_logger, init_logger

PLANET_RADIUS = 6000.0e3
FREQUENCY = 2.0 * np.pi / (86400.0 * 7.5)
UPPER_RADII = (0.2 * PLANET_RADIUS, 0.55 * PLANET_RADIUS, PLANET_RADIUS)

# Three-layer solid-liquid-solid planet. The liquid layer has no shear strength; its shear rheology
# is elastic so that both the classic and the new builders return exactly zero there.
LAYER_KWARGS = dict(
    planet_radius=PLANET_RADIUS,
    forcing_frequency=FREQUENCY,
    density_tuple=(8000.0, 5000.0, 3300.0),
    static_bulk_modulus_tuple=(2.0e11, 1.5e11, 1.0e11),
    static_shear_modulus_tuple=(1.0e11, 0.0, 5.0e10),
    bulk_viscosity_tuple=(1.0e18, 1.0e18, 1.0e18),
    shear_viscosity_tuple=(1.0e20, 1.0e3, 1.0e19),
    layer_type_tuple=("solid", "liquid", "solid"),
    layer_is_static_tuple=(False, True, False),
    layer_is_incompressible_tuple=(False, False, False),
)
SIZE_SPECS = {
    "thickness": dict(thickness_fraction_tuple=(0.2, 0.35, 0.45)),
    "radius": dict(radius_fraction_tuple=(0.2, 0.55, 1.0)),
    "volume": dict(volume_fraction_tuple=(0.2**3, 0.55**3 - 0.2**3, 1.0 - 0.55**3)),
}


def _new_rheologies():
    return dict(shear_rheology_model_tuple=(Maxwell(), Elastic(), Maxwell()),
                bulk_rheology_model_tuple=(Elastic(), Elastic(), Elastic()))


def _classic_rheologies():
    return dict(shear_rheology_model_tuple=(ClassicMaxwell(), ClassicElastic(), ClassicMaxwell()),
                bulk_rheology_model_tuple=(ClassicElastic(), ClassicElastic(), ClassicElastic()))


def _assert_build_data_close(new, classic, modulus_rtol=1.0e-13):
    """Field-by-field comparison of two PlanetBuildData tuples."""
    assert isinstance(new, PlanetBuildData)
    for name in ("radius_array", "density_array", "upper_radius_bylayer_array"):
        np.testing.assert_allclose(getattr(new, name), getattr(classic, name), rtol=1.0e-14, atol=0.0)
    for name in ("complex_bulk_modulus_array", "complex_shear_modulus_array"):
        np.testing.assert_allclose(getattr(new, name), getattr(classic, name), rtol=modulus_rtol, atol=0.0)
    assert new.frequency == classic.frequency
    np.testing.assert_allclose(new.planet_bulk_density, classic.planet_bulk_density, rtol=1.0e-14)
    assert new.layer_types == tuple(classic.layer_types)
    assert new.is_static_bylayer == tuple(classic.is_static_bylayer)
    assert new.is_incompressible_bylayer == tuple(classic.is_incompressible_bylayer)


def _clean_grid_arrays(slices_tuple=(6, 8, 12)):
    """Radially resolved arrays on a solver-ready grid (r = 0, duplicated interfaces, layer tops)."""
    classic = classic_homogeneous(
        slices_tuple=slices_tuple, thickness_fraction_tuple=(0.2, 0.35, 0.45),
        **_classic_rheologies(), **LAYER_KWARGS)
    radius = np.asarray(classic.radius_array)
    layer_index = np.repeat(np.arange(3), slices_tuple)
    per_layer = lambda key: np.asarray(LAYER_KWARGS[key])[layer_index]
    return dict(
        radius_array=radius,
        density_array=per_layer("density_tuple"),
        static_bulk_modulus_array=per_layer("static_bulk_modulus_tuple"),
        static_shear_modulus_array=per_layer("static_shear_modulus_tuple"),
        bulk_viscosity_array=per_layer("bulk_viscosity_tuple"),
        shear_viscosity_array=per_layer("shear_viscosity_tuple"),
    )


def _from_data_kwargs(arrays, **rheologies):
    return dict(
        forcing_frequency=FREQUENCY,
        layer_upper_radius_tuple=UPPER_RADII,
        layer_type_tuple=LAYER_KWARGS["layer_type_tuple"],
        layer_is_static_tuple=LAYER_KWARGS["layer_is_static_tuple"],
        layer_is_incompressible_tuple=LAYER_KWARGS["layer_is_incompressible_tuple"],
        **arrays, **rheologies)


@pytest.fixture
def spdlog_text(tmp_path):
    """Route the C++ logger to a temporary file for the test and hand back a reader for its text."""
    log_path = tmp_path / "tidalpy_x.log"
    init_logger({"console_level": "off", "file_level": "warning", "log_to_file": True,
                 "log_file_path": str(log_path)})

    def read():
        flush_logger()
        return log_path.read_text(encoding="utf-8") if log_path.exists() else ""

    yield read
    init_logger(build_logging_x_config())


# ---------------------------------------------------------------------------------------------------------------------
# Homogeneous-layer builder
# ---------------------------------------------------------------------------------------------------------------------

@pytest.mark.parametrize("slices", (None, (6, 8, 12)), ids=("slice_per_layer", "slices_tuple"))
@pytest.mark.parametrize("size_spec", sorted(SIZE_SPECS))
def test_homogeneous_matches_classic(size_spec, slices):
    """Every output field agrees with the classic builder for each layer-size description."""
    new = build_rs_input_homogeneous_layers(
        slices_tuple=slices, slice_per_layer=7, **SIZE_SPECS[size_spec], **_new_rheologies(), **LAYER_KWARGS)
    classic = classic_homogeneous(
        slices_tuple=slices, slice_per_layer=7, **SIZE_SPECS[size_spec], **_classic_rheologies(), **LAYER_KWARGS)
    _assert_build_data_close(new, classic)


def test_homogeneous_grid_structure():
    """Each layer grid runs from its base to its top so interface radii appear twice."""
    data = build_rs_input_homogeneous_layers(
        slices_tuple=(6, 8, 12), thickness_fraction_tuple=(0.2, 0.35, 0.45), **_new_rheologies(), **LAYER_KWARGS)
    radius = data.radius_array
    assert radius.size == 26
    assert radius[0] == 0.0
    assert radius[-1] == PLANET_RADIUS
    for upper in UPPER_RADII[:-1]:
        assert np.count_nonzero(np.isclose(radius, upper, rtol=1e-12, atol=0.0)) == 2
    np.testing.assert_allclose(data.upper_radius_bylayer_array, UPPER_RADII, rtol=1e-14)
    # Liquid layer: elastic rheology on a zero static shear modulus gives exactly zero.
    assert np.all(data.complex_shear_modulus_array[6:14] == 0.0)
    assert np.all(data.complex_shear_modulus_array[:6].imag > 0.0)
    # Bulk density is the volume-weighted mean of the layer densities.
    r3 = np.asarray((0.0,) + UPPER_RADII) ** 3
    expected = np.sum(np.asarray(LAYER_KWARGS["density_tuple"]) * np.diff(r3)) / PLANET_RADIUS**3
    np.testing.assert_allclose(data.planet_bulk_density, expected, rtol=1e-14)


def test_single_rheology_model_is_broadcast():
    """One model instance (or model name) applies to every layer."""
    per_layer = build_rs_input_homogeneous_layers(
        thickness_fraction_tuple=(0.2, 0.35, 0.45),
        shear_rheology_model_tuple=(Andrade(), Andrade(), Andrade()),
        bulk_rheology_model_tuple=(Elastic(), Elastic(), Elastic()), **LAYER_KWARGS)
    single = build_rs_input_homogeneous_layers(
        thickness_fraction_tuple=(0.2, 0.35, 0.45),
        shear_rheology_model_tuple=Andrade(), bulk_rheology_model_tuple=Elastic(), **LAYER_KWARGS)
    by_name = build_rs_input_homogeneous_layers(
        thickness_fraction_tuple=(0.2, 0.35, 0.45),
        shear_rheology_model_tuple="andrade", bulk_rheology_model_tuple=("elastic", Elastic(), "elastic"),
        **LAYER_KWARGS)
    for other in (single, by_name):
        np.testing.assert_array_equal(other.complex_shear_modulus_array, per_layer.complex_shear_modulus_array)
        np.testing.assert_array_equal(other.complex_bulk_modulus_array, per_layer.complex_bulk_modulus_array)


def test_rejects_classic_rheology_models():
    """Classic TidalPy.rheology models raise a TypeError that points at rheology_x."""
    with pytest.raises(TypeError, match="rheology_x"):
        build_rs_input_homogeneous_layers(
            thickness_fraction_tuple=(0.2, 0.35, 0.45),
            shear_rheology_model_tuple=ClassicMaxwell(), bulk_rheology_model_tuple=Elastic(), **LAYER_KWARGS)
    with pytest.raises(TypeError, match="entry 1"):
        build_rs_input_homogeneous_layers(
            thickness_fraction_tuple=(0.2, 0.35, 0.45),
            shear_rheology_model_tuple=(Maxwell(), ClassicMaxwell(), Maxwell()),
            bulk_rheology_model_tuple=Elastic(), **LAYER_KWARGS)


@pytest.mark.parametrize("override, match", (
    (dict(static_bulk_modulus_tuple=(2.0e11, 1.5e11)), "static_bulk_modulus_tuple"),
    (dict(shear_rheology_model_tuple=(Maxwell(), Maxwell())), "one rheology model per layer"),
    (dict(thickness_fraction_tuple=(0.2, 0.35, 0.40)), "sum to 1"),
    (dict(thickness_fraction_tuple=None, radius_fraction_tuple=(0.2, 0.15, 1.0)), "increase"),
    (dict(thickness_fraction_tuple=None, volume_fraction_tuple=(0.5, -0.1, 0.6)), "positive"),
    (dict(slices_tuple=(6, 4, 12)), "at least 5"),
    (dict(radius_fraction_tuple=(0.2, 0.55, 1.0)), "exactly one"),
    (dict(thickness_fraction_tuple=None), "exactly one"),
))
def test_homogeneous_invalid_inputs_raise_value_error(override, match):
    """Bad sizes, fractions, or slice counts raise ValueError with a pointed message."""
    kwargs = dict(thickness_fraction_tuple=(0.2, 0.35, 0.45), **_new_rheologies(), **LAYER_KWARGS)
    kwargs.update(override)
    with pytest.raises(ValueError, match=match):
        build_rs_input_homogeneous_layers(**kwargs)


# ---------------------------------------------------------------------------------------------------------------------
# From-data builder
# ---------------------------------------------------------------------------------------------------------------------

def test_from_data_matches_classic_on_clean_grid():
    """On a solver-ready grid the new builder reproduces the classic one field by field."""
    arrays = _clean_grid_arrays()
    new = build_rs_input_from_data(**_from_data_kwargs(arrays, **_new_rheologies()), warnings=False)
    classic = classic_from_data(**_from_data_kwargs(arrays, **_classic_rheologies()), warnings=False)
    _assert_build_data_close(new, classic)
    np.testing.assert_array_equal(new.radius_array, arrays["radius_array"])


def test_from_data_repairs_grid_and_warns(spdlog_text):
    """A grid missing r = 0, a whole interface, and an interface duplicate is repaired with a warning each."""
    arrays = _clean_grid_arrays(slices_tuple=(7, 8, 12))
    clean_radius = arrays["radius_array"]
    # Clean layout: layer 0 = slices 0-6, layer 1 = 7-14, layer 2 = 15-26. Drop r = 0, both copies of the
    # first interface (layer 0 loses its top, layer 1 its base), and layer 2's copy of the second
    # interface. A singly listed interface slice is taken as the top of the lower layer (the classic
    # builder's convention), so the upper layer's base is re-inserted with the next slice's properties.
    keep = np.ones(clean_radius.size, dtype=bool)
    keep[[0, 6, 7, 15]] = False
    broken = {key: value[keep] for key, value in arrays.items()}

    data = build_rs_input_from_data(**_from_data_kwargs(broken, **_new_rheologies()), warnings=True)
    radius = data.radius_array
    assert radius[0] == 0.0
    assert radius[-1] == PLANET_RADIUS
    assert radius.size == clean_radius.size
    np.testing.assert_allclose(radius, clean_radius, rtol=1e-14, atol=0.0)
    for upper in UPPER_RADII[:-1]:
        assert np.count_nonzero(np.isclose(radius, upper, rtol=1e-12, atol=0.0)) == 2
    # Inserted slices copy their neighbour's properties, so the repaired arrays equal the clean ones.
    np.testing.assert_array_equal(data.density_array, arrays["density_array"])
    reference = build_rs_input_from_data(**_from_data_kwargs(arrays, **_new_rheologies()), warnings=False)
    np.testing.assert_array_equal(data.complex_shear_modulus_array, reference.complex_shear_modulus_array)
    np.testing.assert_array_equal(data.complex_bulk_modulus_array, reference.complex_bulk_modulus_array)
    # The bulk density is the piecewise-constant shell mass over the planet volume (the classic builder
    # mis-attributes the shell below an inserted layer top, so it is not the reference here).
    shell_volume = (4.0 / 3.0) * np.pi * np.diff(radius**3)
    expected = np.sum(shell_volume * data.density_array[1:]) / ((4.0 / 3.0) * np.pi * PLANET_RADIUS**3)
    np.testing.assert_allclose(data.planet_bulk_density, expected, rtol=1e-12)
    np.testing.assert_allclose(data.planet_bulk_density, reference.planet_bulk_density, rtol=1e-12)

    # NOTE (0.9.0): the classic builder repairs the arrays the same way.
    classic = classic_from_data(**_from_data_kwargs(broken, **_classic_rheologies()), warnings=False)
    np.testing.assert_allclose(classic.radius_array, radius, rtol=1e-14, atol=0.0)
    np.testing.assert_array_equal(classic.density_array, data.density_array)

    text = spdlog_text()
    assert text.count("build_rs_input_from_data") == 4
    assert "start at zero" in text
    assert text.count("appear twice") == 2
    assert "does not have its upper radius" in text


def test_from_data_warnings_can_be_silenced(spdlog_text):
    """`warnings=False` repairs the grid without logging."""
    arrays = _clean_grid_arrays()
    broken = {key: value[1:] for key, value in arrays.items()}  # drop r = 0
    data = build_rs_input_from_data(**_from_data_kwargs(broken, **_new_rheologies()), warnings=False)
    assert data.radius_array[0] == 0.0
    assert "build_rs_input_from_data" not in spdlog_text()


def test_from_data_accepts_lists_and_single_model():
    """Plain lists are accepted for the arrays and a single model is applied to every layer."""
    arrays = {key: list(value) for key, value in _clean_grid_arrays().items()}
    data = build_rs_input_from_data(
        **_from_data_kwargs(arrays, shear_rheology_model_tuple=Maxwell(), bulk_rheology_model_tuple="elastic"),
        warnings=False)
    reference = build_rs_input_from_data(
        **_from_data_kwargs(_clean_grid_arrays(), shear_rheology_model_tuple=(Maxwell(),) * 3,
                            bulk_rheology_model_tuple=(Elastic(),) * 3), warnings=False)
    np.testing.assert_array_equal(data.complex_shear_modulus_array, reference.complex_shear_modulus_array)
    np.testing.assert_array_equal(data.complex_bulk_modulus_array, reference.complex_bulk_modulus_array)


@pytest.mark.parametrize("mutate, match", (
    (lambda a, k: a.__setitem__("density_array", a["density_array"][:-1]), "same length"),
    (lambda a, k: a["radius_array"].__setitem__(3, a["radius_array"][5]), "ascending"),
    (lambda a, k: k.__setitem__("layer_upper_radius_tuple", (UPPER_RADII[0], UPPER_RADII[1], 0.9 * PLANET_RADIUS)),
     "planet radius"),
    (lambda a, k: k.__setitem__("layer_upper_radius_tuple", (UPPER_RADII[1], UPPER_RADII[0], PLANET_RADIUS)),
     "increase"),
    (lambda a, k: k.__setitem__("layer_type_tuple", ("solid", "liquid")), "layer_type_tuple"),
))
def test_from_data_invalid_inputs_raise_value_error(mutate, match):
    """Inconsistent arrays or layer boundaries raise ValueError with a pointed message."""
    arrays = {key: np.array(value) for key, value in _clean_grid_arrays().items()}
    kwargs = _from_data_kwargs(arrays, **_new_rheologies())
    mutate(kwargs, kwargs)
    with pytest.raises(ValueError, match=match):
        build_rs_input_from_data(**kwargs, warnings=False)


def test_from_data_rejects_classic_rheology_models():
    arrays = _clean_grid_arrays()
    with pytest.raises(TypeError, match="rheology_x"):
        build_rs_input_from_data(
            **_from_data_kwargs(arrays, shear_rheology_model_tuple=ClassicMaxwell(),
                                bulk_rheology_model_tuple=Elastic()), warnings=False)


# ---------------------------------------------------------------------------------------------------------------------
# Output feeds the solvers
# ---------------------------------------------------------------------------------------------------------------------

def test_output_feeds_both_solvers():
    """The named tuple unpacks positionally into either radial solver."""
    data = build_rs_input_homogeneous_layers(
        slices_tuple=(8, 8, 12), thickness_fraction_tuple=(0.2, 0.35, 0.45), **_new_rheologies(), **LAYER_KWARGS)
    solver_kwargs = dict(degree_l=2, solve_for=("tidal",), integration_rtol=1.0e-8, integration_atol=1.0e-10)
    new = radial_solver(*data, **solver_kwargs)
    assert new.success, new.message
    k_new = complex(np.atleast_1d(new.k)[0])
    assert np.isfinite(k_new) and 0.0 < k_new.real < 1.5
    old = radial_solver_old(*data, **solver_kwargs)   # NOTE (0.9.0): classic solver reference.
    assert old.success, old.message
    np.testing.assert_allclose(np.atleast_1d(old.k), np.atleast_1d(new.k), rtol=1e-6)
