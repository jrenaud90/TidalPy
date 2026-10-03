"""The radial solver input builders against the frozen classic builders, plus their validation and grid repair."""
from pathlib import Path

import numpy as np
import pytest

from TidalPy.RadialSolver import (
    PlanetBuildData,
    build_rs_input_from_data,
    build_rs_input_homogeneous_layers,
    radial_solver,
)
from TidalPy.Rheology import Andrade, Elastic, Maxwell

# Frozen classic builder outputs keyed '<case>__<PlanetBuildData field>', plus the classic solver's k.
FROZEN_PATH = Path(__file__).parent / "frozen" / "test_l_build_inputs.npz"
with np.load(FROZEN_PATH, allow_pickle=False) as frozen_file:
    CLASSIC_REFERENCE = {key: frozen_file[key] for key in frozen_file.files}

PLANET_RADIUS = 6000.0e3
FREQUENCY = 2.0 * np.pi / (86400.0 * 7.5)
UPPER_RADII = (0.2 * PLANET_RADIUS, 0.55 * PLANET_RADIUS, PLANET_RADIUS)
THICKNESS_FRACTIONS = (0.2, 0.35, 0.45)

# Solid-liquid-solid planet. The liquid's shear rheology is elastic so both builders return exactly zero there.
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
    "thickness": dict(thickness_fraction_tuple=THICKNESS_FRACTIONS),
    "radius": dict(radius_fraction_tuple=(0.2, 0.55, 1.0)),
    "volume": dict(volume_fraction_tuple=(0.2**3, 0.55**3 - 0.2**3, 1.0 - 0.55**3)),
}


def _new_rheologies():
    return dict(shear_rheology_model_tuple=(Maxwell(), Elastic(), Maxwell()),
                bulk_rheology_model_tuple=(Elastic(), Elastic(), Elastic()))


def _build_homogeneous(**kwargs):
    """The homogeneous builder on the three-layer planet with thickness fractions and the Maxwell rheologies."""
    return build_rs_input_homogeneous_layers(
        **(dict(thickness_fraction_tuple=THICKNESS_FRACTIONS) | _new_rheologies() | kwargs),
        **LAYER_KWARGS)


def _frozen_classic_build_data(case):
    """The classic builder's frozen output for `case` as a PlanetBuildData tuple."""
    def field(name):
        return CLASSIC_REFERENCE[f"{case}__{name}"]

    return PlanetBuildData(
        radius_array=field("radius_array"),
        density_array=field("density_array"),
        complex_bulk_modulus_array=field("complex_bulk_modulus_array"),
        complex_shear_modulus_array=field("complex_shear_modulus_array"),
        frequency=float(field("frequency")),
        planet_bulk_density=float(field("planet_bulk_density")),
        layer_types=tuple(str(layer_type) for layer_type in field("layer_types")),
        is_static_bylayer=tuple(bool(flag) for flag in field("is_static_bylayer")),
        is_incompressible_bylayer=tuple(bool(flag) for flag in field("is_incompressible_bylayer")),
        upper_radius_bylayer_array=field("upper_radius_bylayer_array"),
    )


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
    grid = _build_homogeneous(slices_tuple=slices_tuple)
    layer_index = np.repeat(np.arange(3), slices_tuple)

    def per_layer(key):
        return np.asarray(LAYER_KWARGS[key])[layer_index]

    return dict(
        radius_array=np.asarray(grid.radius_array),
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


# Homogeneous-layer builder

@pytest.mark.parametrize("slices", (None, (6, 8, 12)), ids=("slice_per_layer", "slices_tuple"))
@pytest.mark.parametrize("size_spec", sorted(SIZE_SPECS))
def test_homogeneous_matches_classic(size_spec, slices):
    """Every output field agrees with the classic builder for each layer-size description."""
    new = build_rs_input_homogeneous_layers(
        slices_tuple=slices,
        slice_per_layer=7,
        **SIZE_SPECS[size_spec],
        **_new_rheologies(),
        **LAYER_KWARGS)
    slices_id = "slice_per_layer" if slices is None else "slices_tuple"
    classic = _frozen_classic_build_data(f"homogeneous__{size_spec}__{slices_id}")
    _assert_build_data_close(new, classic)


def test_homogeneous_grid_structure():
    """Interface radii appear twice, the liquid has zero shear, and the bulk density is the volume-weighted mean."""
    data = _build_homogeneous(slices_tuple=(6, 8, 12))
    radius = data.radius_array
    assert radius.size == 26
    assert radius[0] == 0.0
    assert radius[-1] == PLANET_RADIUS
    for upper in UPPER_RADII[:-1]:
        assert np.count_nonzero(np.isclose(radius, upper, rtol=1e-12, atol=0.0)) == 2
    np.testing.assert_allclose(data.upper_radius_bylayer_array, UPPER_RADII, rtol=1e-14)
    assert np.all(data.complex_shear_modulus_array[6:14] == 0.0)
    assert np.all(data.complex_shear_modulus_array[:6].imag > 0.0)
    radius_cubed = np.asarray((0.0,) + UPPER_RADII) ** 3
    expected = np.sum(np.asarray(LAYER_KWARGS["density_tuple"]) * np.diff(radius_cubed)) / PLANET_RADIUS**3
    np.testing.assert_allclose(data.planet_bulk_density, expected, rtol=1e-14)


def test_single_rheology_model_is_broadcast():
    """One model instance (or model name) applies to every layer."""
    per_layer = _build_homogeneous(
        shear_rheology_model_tuple=(Andrade(), Andrade(), Andrade()),
        bulk_rheology_model_tuple=(Elastic(), Elastic(), Elastic()))
    single = _build_homogeneous(shear_rheology_model_tuple=Andrade(), bulk_rheology_model_tuple=Elastic())
    by_name = _build_homogeneous(
        shear_rheology_model_tuple="andrade", bulk_rheology_model_tuple=("elastic", Elastic(), "elastic"))
    for other in (single, by_name):
        np.testing.assert_array_equal(other.complex_shear_modulus_array, per_layer.complex_shear_modulus_array)
        np.testing.assert_array_equal(other.complex_bulk_modulus_array, per_layer.complex_bulk_modulus_array)


@pytest.mark.parametrize("shear_models, match", (
    (object(), r"`shear_rheology_model_tuple` must be"),
    ((Andrade(), object(), Andrade()), r"`shear_rheology_model_tuple` entry 1 must be"),
))
def test_rejects_objects_that_are_not_rheology_models(shear_models, match):
    """Anything but a Rheology model or a model name raises a TypeError naming the argument and tuple entry."""
    with pytest.raises(TypeError, match=match):
        _build_homogeneous(shear_rheology_model_tuple=shear_models, bulk_rheology_model_tuple=Elastic())


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
    kwargs = dict(thickness_fraction_tuple=THICKNESS_FRACTIONS, **_new_rheologies(), **LAYER_KWARGS)
    kwargs.update(override)
    with pytest.raises(ValueError, match=match):
        build_rs_input_homogeneous_layers(**kwargs)


# From-data builder

def test_from_data_matches_classic_on_clean_grid():
    """On a solver-ready grid the new builder reproduces the classic one field by field."""
    arrays = _clean_grid_arrays()
    new = build_rs_input_from_data(**_from_data_kwargs(arrays, **_new_rheologies()), warnings=False)
    _assert_build_data_close(new, _frozen_classic_build_data("from_data_clean"))
    np.testing.assert_array_equal(new.radius_array, arrays["radius_array"])


def test_from_data_repairs_grid_and_warns(spdlog_text):
    """A grid missing r = 0, a whole interface, and an interface duplicate is repaired with a warning each."""
    arrays = _clean_grid_arrays(slices_tuple=(7, 8, 12))
    clean_radius = arrays["radius_array"]
    # Layers span slices 0-6, 7-14, and 15-26. Drop r = 0, both copies of the first interface, and layer 2's copy of
    # the second. A singly listed interface is the lower layer's top, so the upper base is re-inserted from above.
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
    # Inserted slices copy their neighbor's properties, so the repaired arrays equal the clean ones.
    np.testing.assert_array_equal(data.density_array, arrays["density_array"])
    reference = build_rs_input_from_data(**_from_data_kwargs(arrays, **_new_rheologies()), warnings=False)
    np.testing.assert_array_equal(data.complex_shear_modulus_array, reference.complex_shear_modulus_array)
    np.testing.assert_array_equal(data.complex_bulk_modulus_array, reference.complex_bulk_modulus_array)
    # The classic builder mis-attributes the shell below an inserted layer top, so the expected bulk density is
    # computed directly rather than taken from it.
    shell_volume = (4.0 / 3.0) * np.pi * np.diff(radius**3)
    expected = np.sum(shell_volume * data.density_array[1:]) / ((4.0 / 3.0) * np.pi * PLANET_RADIUS**3)
    np.testing.assert_allclose(data.planet_bulk_density, expected, rtol=1e-12)
    np.testing.assert_allclose(data.planet_bulk_density, reference.planet_bulk_density, rtol=1e-12)

    np.testing.assert_allclose(CLASSIC_REFERENCE["from_data_repaired__radius_array"], radius, rtol=1e-14, atol=0.0)
    np.testing.assert_array_equal(CLASSIC_REFERENCE["from_data_repaired__density_array"], data.density_array)

    text = spdlog_text()
    assert text.count("build_rs_input_from_data") == 4
    assert "start at zero" in text
    assert text.count("appear twice") == 2
    assert "does not have its upper radius" in text


def test_from_data_warnings_can_be_silenced(spdlog_text):
    """`warnings=False` repairs the grid without logging."""
    broken = {key: value[1:] for key, value in _clean_grid_arrays().items()}
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
    (lambda kwargs: kwargs.__setitem__("density_array", kwargs["density_array"][:-1]), "same length"),
    (lambda kwargs: kwargs["radius_array"].__setitem__(3, kwargs["radius_array"][5]), "ascending"),
    (lambda kwargs: kwargs.__setitem__("layer_upper_radius_tuple", UPPER_RADII[:2] + (0.9 * PLANET_RADIUS,)),
     "planet radius"),
    (lambda kwargs: kwargs.__setitem__("layer_upper_radius_tuple", (UPPER_RADII[1], UPPER_RADII[0], PLANET_RADIUS)),
     "increase"),
    (lambda kwargs: kwargs.__setitem__("layer_type_tuple", ("solid", "liquid")), "layer_type_tuple"),
), ids=("density_length", "radius_order", "top_below_surface", "upper_radius_order", "layer_count"))
def test_from_data_invalid_inputs_raise_value_error(mutate, match):
    """Inconsistent arrays or layer boundaries raise ValueError with a pointed message."""
    arrays = {key: np.array(value) for key, value in _clean_grid_arrays().items()}
    kwargs = _from_data_kwargs(arrays, **_new_rheologies())
    mutate(kwargs)
    with pytest.raises(ValueError, match=match):
        build_rs_input_from_data(**kwargs, warnings=False)


def test_output_feeds_the_solver():
    """The named tuple unpacks positionally into the radial solver, which reproduces the classic solver's k."""
    data = _build_homogeneous(slices_tuple=(8, 8, 12))
    new = radial_solver(
        *data,
        degree_l=2,
        solve_for=("tidal",),
        integration_rtol=1.0e-8,
        integration_atol=1.0e-10)
    assert new.success, new.message
    k_new = complex(np.atleast_1d(new.k)[0])
    assert np.isfinite(k_new) and 0.0 < k_new.real < 1.5
    np.testing.assert_allclose(CLASSIC_REFERENCE["classic_solver_k"], np.atleast_1d(new.k), rtol=1e-6)
