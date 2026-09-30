"""Material EOS models (constant, Birch-Murnaghan, Vinet, interpolated): density, inversion, thermal terms, I/O."""

import math
import struct

import numpy as np
import pytest

from TidalPy.Material.eos import material_eos
from TidalPy.Material.eos import make_material_eos
from TidalPy.Utilities.classes.classes import PhysicsBase, TidalPyBaseClass

_RHO0 = 3500.0      # [kg/m^3]
_K0 = 1.30e11       # [Pa]
_K0P = 4.5          # [dimensionless]
_ALPHA = 3.0e-5     # [1/K]
_T_REF = 300.0      # [K]

_LAWS = [("BirchMurnaghanEOS", "birch_murnaghan_pressure"), ("VinetEOS", "vinet_pressure")]
_MODEL_NAMES = ["ConstantDensityEOS", "BirchMurnaghanEOS", "VinetEOS", "InterpolatedEOS"]


def _binary_copy(model, blank, tmp_path):
    """Save a model to a binary file and load it into a blank model."""
    path = str(tmp_path / "eos.tpyb")
    model.save_binary(path)
    blank.load_binary(path)
    return blank


def _thermal_model(cls_name):
    if cls_name == "InterpolatedEOS":
        return material_eos.InterpolatedEOS([0.0, 1.0e6, 2.0e6], [5000.0, 4000.0, 3000.0], thermal_expansion=_ALPHA)
    if cls_name == "ConstantDensityEOS":
        return material_eos.ConstantDensityEOS(_RHO0, thermal_expansion=_ALPHA)
    return getattr(material_eos, cls_name)(_RHO0, _K0, _K0P, thermal_expansion=_ALPHA)


def test_constant_density():
    """Constant density ignores pressure and radius."""
    eos = material_eos.ConstantDensityEOS(reference_density=_RHO0)
    assert eos.reference_density == pytest.approx(_RHO0)
    assert eos.calc_density(0.0) == pytest.approx(_RHO0)
    assert eos.calc_density(5.0e10) == pytest.approx(_RHO0)
    assert eos.calc_density(0.0, 0.0, 1.0e6) == pytest.approx(_RHO0)


# =====================================================================================================================
# Analytic pressure laws and inversion
# =====================================================================================================================
@pytest.mark.parametrize("cls_name,law", _LAWS)
def test_reference_state_has_zero_pressure(cls_name, law):
    """The law is zero at eta = 1 and the density at zero pressure is the reference density."""
    assert getattr(material_eos, law)(1.0, _K0, _K0P) == pytest.approx(0.0, abs=1e-3)
    eos = getattr(material_eos, cls_name)(_RHO0, _K0, _K0P)
    assert eos.calc_density(0.0) == pytest.approx(_RHO0, rel=1e-9)


@pytest.mark.parametrize("cls_name,law", _LAWS)
@pytest.mark.parametrize("pressure", [1.0e9, 1.0e10, 5.0e10, 1.5e11])
def test_inversion_roundtrip(cls_name, law, pressure):
    """The inverted density pushed back through the forward law recovers the pressure."""
    rho = getattr(material_eos, cls_name)(_RHO0, _K0, _K0P).calc_density(pressure)
    assert rho > _RHO0
    assert getattr(material_eos, law)(rho / _RHO0, _K0, _K0P) == pytest.approx(pressure, rel=1e-6)


@pytest.mark.parametrize("cls_name", ["BirchMurnaghanEOS", "VinetEOS"])
def test_density_monotonic_in_pressure(cls_name):
    eos = getattr(material_eos, cls_name)(_RHO0, _K0, _K0P)
    rhos = [eos.calc_density(pressure) for pressure in (0.0, 1e10, 5e10, 1e11)]
    assert all(later > earlier for earlier, later in zip(rhos, rhos[1:]))


@pytest.mark.parametrize("cls_name,law", _LAWS)
@pytest.mark.parametrize("bulk_modulus_derivative", [3.2, 4.0, 5.5])
@pytest.mark.parametrize("pressure", [-5.0e9, -1.0e6, 1.0e3, 1.0e9, 1.4e11, 2.5e11])
def test_inversion_roundtrip_is_tight_in_compression_and_tension(cls_name, law, bulk_modulus_derivative, pressure):
    """Inside the monotonic range the inversion recovers the pressure to its tolerance."""
    eos = getattr(material_eos, cls_name)(_RHO0, _K0, bulk_modulus_derivative)
    eta = eos.calc_density(pressure) / _RHO0
    assert (eta > 1.0) == (pressure > 0.0)
    recovered = getattr(material_eos, law)(eta, _K0, bulk_modulus_derivative)
    # K is the slope, so a compression good to 1e-13 gives a pressure good to about K 1e-13.
    assert recovered == pytest.approx(pressure, rel=1.0e-10, abs=_K0 * 1.0e-11)


@pytest.mark.parametrize("cls_name,law", _LAWS)
@pytest.mark.parametrize("bulk_modulus_derivative", [3.2, 4.0, 5.5])
def test_density_is_continuous_where_the_law_turns_over_in_tension(cls_name, law, bulk_modulus_derivative):
    """Past the law's minimum pressure the density holds the turning point, with no jump."""
    # The structure solve probes deep tension while its central pressure is a guess; a jump stalls the stepper.
    law_func = getattr(material_eos, law)
    eos = getattr(material_eos, cls_name)(_RHO0, _K0, bulk_modulus_derivative)
    floor_density = eos.calc_density(-10.0 * _K0)
    assert eos.calc_density(-100.0 * _K0) == floor_density
    floor_eta = floor_density / _RHO0
    assert 0.0 < floor_eta < 1.0
    floor_pressure = law_func(floor_eta, _K0, bulk_modulus_derivative)
    assert law_func(0.99 * floor_eta, _K0, bulk_modulus_derivative) > floor_pressure
    # The law is flat at its minimum, so 1e-9 K0 in pressure is about its square root in compression.
    just_inside = eos.calc_density(floor_pressure + 1.0e-9 * _K0)
    assert just_inside == pytest.approx(floor_density, rel=1.0e-3)
    assert just_inside >= floor_density


def test_birch_murnaghan_holds_its_maximum_pressure_when_it_turns_over_in_compression():
    """With K0' < 4 the law turns over at large compression and the density holds; with K0' >= 4 it keeps rising."""
    eos = material_eos.BirchMurnaghanEOS(_RHO0, _K0, 3.2)
    ceiling_density = eos.calc_density(1.0e15)
    assert eos.calc_density(1.0e16) == ceiling_density
    ceiling_eta = ceiling_density / _RHO0
    ceiling_pressure = material_eos.birch_murnaghan_pressure(ceiling_eta, _K0, 3.2)
    assert material_eos.birch_murnaghan_pressure(1.01 * ceiling_eta, _K0, 3.2) < ceiling_pressure
    unbounded = material_eos.BirchMurnaghanEOS(_RHO0, _K0, 4.5)
    assert unbounded.calc_density(1.0e16) > unbounded.calc_density(1.0e15)


@pytest.mark.parametrize("cls_name", ["BirchMurnaghanEOS", "VinetEOS"])
def test_a_reloaded_model_inverts_like_the_one_that_was_saved(cls_name, tmp_path):
    """The monotonic range is derived from K0 and K0', so a binary load must rebuild it."""
    eos = getattr(material_eos, cls_name)(_RHO0, 2.2e11, 3.4)
    reloaded = _binary_copy(eos, getattr(material_eos, cls_name)(_RHO0, _K0, _K0P), tmp_path)
    for pressure in (-1.0e13, -1.0e9, 5.0e10, 1.0e15):
        assert reloaded.calc_density(pressure) == eos.calc_density(pressure)


def test_bm_and_vinet_agree_at_small_compression():
    birch_murnaghan = material_eos.BirchMurnaghanEOS(_RHO0, _K0, _K0P)
    vinet = material_eos.VinetEOS(_RHO0, _K0, _K0P)
    assert birch_murnaghan.calc_density(5.0e9) == pytest.approx(vinet.calc_density(5.0e9), rel=0.02)


# =====================================================================================================================
# InterpolatedEOS
# =====================================================================================================================
@pytest.mark.parametrize("radius,expected", [(0.0, 5000.0), (2.0e6, 3000.0), (0.5e6, 4500.0), (-1.0e6, 5000.0),
                                             (9.0e6, 3000.0)],
                         ids=["first_node", "last_node", "midpoint", "clamp_below", "clamp_above"])
def test_interpolated_density(radius, expected):
    eos = material_eos.InterpolatedEOS([0.0, 1.0e6, 2.0e6], [5000.0, 4000.0, 3000.0])
    assert eos.num_points == 3
    assert eos.calc_density(0.0, 0.0, radius) == pytest.approx(expected)


def test_interpolated_reads_match_numpy_on_an_uneven_table(tmp_path):
    """Every table of a PREM-like table (uneven, with a discontinuity) reads as numpy.interp, before and after I/O."""
    radius = np.concatenate([np.linspace(0.0, 3.0e6, 7), [3.0e6], np.linspace(3.2e6, 6.0e6, 9),
                             6.0e6 + np.cumsum(np.geomspace(1.0e5, 2.0e3, 40))])
    radius = np.sort(radius)
    rows = np.arange(radius.size, dtype=float)
    tables = {"density": 1.0e4 - 10.0 * rows, "shear_modulus": 5.0e10 + 1.0e8 * rows ** 1.5,
              "bulk_modulus": 1.0e11 + 3.0e8 * rows, "shear_viscosity": 10.0 ** (18.0 + 0.05 * rows),
              "bulk_viscosity": 10.0 ** (20.0 - 0.02 * rows)}
    # Distinct values across the discontinuity, so reading the wrong side of it would show.
    discontinuity = int(np.flatnonzero(np.diff(radius) == 0.0)[0])
    for values in tables.values():
        values[discontinuity + 1:] *= 0.7
    eos = material_eos.InterpolatedEOS(radius, **tables)

    queries = np.concatenate([radius, radius[:-1] + 0.5 * np.diff(radius), np.nextafter(radius, -np.inf),
                              np.linspace(-1.0e6, radius[-1] + 1.0e6, 3001)])
    readers = {"density": lambda model, r: model.calc_density(0.0, None, r),
               "shear_modulus": lambda model, r: model.get_tabulated_shear_modulus(r),
               "bulk_modulus": lambda model, r: model.get_tabulated_bulk_modulus(r),
               "shear_viscosity": lambda model, r: model.get_tabulated_shear_viscosity(r),
               "bulk_viscosity": lambda model, r: model.get_tabulated_bulk_viscosity(r)}

    def check(model):
        for name, read in readers.items():
            expected = np.interp(queries, radius, tables[name])
            got = np.array([read(model, r) for r in queries])
            np.testing.assert_allclose(
                got,
                expected,
                rtol=1.0e-14,
                atol=0.0,
                err_msg=name)

    check(eos)
    # A binary load rebuilds the search seeds.
    check(_binary_copy(eos, make_material_eos(eos.model_name), tmp_path))


def _saved_interpolated_bytes(tmp_path):
    """A three-point interpolated EOS saved to a file, returned as its bytes."""
    path = tmp_path / "eos.tpyb"
    material_eos.InterpolatedEOS([0.0, 1.0e6, 2.0e6], [5000.0, 4000.0, 3000.0]).save_binary(str(path))
    return path.read_bytes()


def _load_refused(data, tmp_path, match):
    """Loading the bytes raises IOError matching ``match`` and leaves the target model as it was."""
    path = tmp_path / "bad.tpyb"
    path.write_bytes(bytes(data))
    blank = make_material_eos("interpolate")
    config_before = blank.get_config_dict()
    with pytest.raises(IOError, match=match):
        blank.load_binary(str(path))
    assert blank.get_config_dict() == config_before


def test_an_interpolated_table_out_of_order_is_refused_on_load(tmp_path):
    """A radius table out of order fails the load, as it fails the constructor."""
    data = _saved_interpolated_bytes(tmp_path)
    middle_radius = struct.pack("<d", 1.0e6)
    assert data.count(middle_radius) == 1
    _load_refused(data.replace(middle_radius, struct.pack("<d", 3.0e6)), tmp_path, "ascending")


def test_an_interpolated_record_of_the_wrong_size_is_refused_on_load(tmp_path):
    """A header payload size that disagrees with the tables fails the load, as for the other EOS models."""
    data = bytearray(_saved_interpolated_bytes(tmp_path))
    # The payload size is the header's last field: 4 magic bytes, 4 version and byte-order bytes, a 4-byte class id.
    payload_size = struct.unpack_from("<Q", data, 12)[0]
    struct.pack_into("<Q", data, 12, payload_size + 8)
    _load_refused(data + bytes(8), tmp_path, "interpolated EOS record holds")


def test_config_dict_interpolated_roundtrip():
    """An interpolated EOS emits its tables (optional ones only when supplied) and rebuilds through the factory."""
    radii = [0.0, 1.0e6, 2.0e6]
    densities = [5000.0, 4000.0, 3000.0]
    eos = material_eos.InterpolatedEOS(radii, densities)
    cfg = eos.get_config_dict()
    assert cfg["model"] == eos.model_name
    assert cfg["radius_m"] == pytest.approx(radii)
    assert cfg["density_kg_m3"] == pytest.approx(densities)
    for optional_key in ("shear_modulus_pa", "bulk_modulus_pa", "shear_viscosity_pas", "bulk_viscosity_pas"):
        assert optional_key not in cfg
    rebuilt = make_material_eos(cfg["model"], {key: value for key, value in cfg.items() if key != "model"})
    assert rebuilt.num_points == 3
    assert rebuilt.get_config_dict() == cfg


# =====================================================================================================================
# Factory, config dict, and binary round trip
# =====================================================================================================================
@pytest.mark.parametrize("name,cls_attr", [
    ("constant", "ConstantDensityEOS"),
    ("uniform", "ConstantDensityEOS"),
    ("bm", "BirchMurnaghanEOS"),
    ("birch_murnaghan", "BirchMurnaghanEOS"),
    ("vinet", "VinetEOS"),
    ("interp", "InterpolatedEOS"),
])
def test_factory_aliases(name, cls_attr):
    cfg = {"reference_density_kg_m3": _RHO0, "reference_bulk_modulus_pa": _K0,
           "bulk_modulus_derivative": _K0P, "radius_m": [0.0, 1.0e6],
           "density_kg_m3": [5000.0, 4000.0]}
    assert isinstance(make_material_eos(name, cfg), getattr(material_eos, cls_attr))


def test_factory_returns_usable_model():
    eos = make_material_eos("bm", {"reference_density_kg_m3": _RHO0,
                                   "reference_bulk_modulus_pa": _K0,
                                   "bulk_modulus_derivative": _K0P})
    assert eos.calc_density(0.0) == pytest.approx(_RHO0, rel=1e-9)
    assert eos.reference_bulk_modulus == pytest.approx(_K0)


def test_config_dict_bm():
    cfg = material_eos.BirchMurnaghanEOS(_RHO0, _K0, _K0P).get_config_dict()
    assert cfg["model"] == "birch_murnaghan"
    assert cfg["reference_density_kg_m3"] == pytest.approx(_RHO0)
    assert cfg["reference_bulk_modulus_pa"] == pytest.approx(_K0)
    assert cfg["bulk_modulus_derivative"] == pytest.approx(_K0P)
    assert cfg["invert_rtol"] > 0.0
    assert cfg["invert_max_iters"] > 0


@pytest.mark.parametrize("cls_name", ["BirchMurnaghanEOS", "VinetEOS"])
def test_invert_settings_configurable(cls_name):
    """Inversion settings have positive defaults, accept overrides, and a tighter tolerance barely moves the result."""
    model_class = getattr(material_eos, cls_name)
    default = model_class(_RHO0, _K0, _K0P)
    assert default.invert_rtol > 0.0
    assert default.invert_max_iters > 0
    tuned = model_class(
        _RHO0,
        _K0,
        _K0P,
        invert_rtol=1.0e-8,
        invert_max_iters=80)
    assert tuned.invert_rtol == pytest.approx(1.0e-8)
    assert tuned.invert_max_iters == 80
    assert tuned.calc_density(5.0e10) == pytest.approx(default.calc_density(5.0e10), rel=1e-6)


def test_invert_settings_via_factory_and_binary(tmp_path):
    eos = make_material_eos("bm", {"reference_density_kg_m3": _RHO0,
                                   "reference_bulk_modulus_pa": _K0,
                                   "bulk_modulus_derivative": _K0P,
                                   "invert_rtol": 1.0e-9,
                                   "invert_max_iters": 75})
    for model in (eos, _binary_copy(eos, make_material_eos("bm"), tmp_path)):
        assert model.invert_rtol == pytest.approx(1.0e-9)
        assert model.invert_max_iters == 75


@pytest.mark.parametrize("factory", [
    lambda: material_eos.ConstantDensityEOS(_RHO0),
    lambda: material_eos.BirchMurnaghanEOS(_RHO0, _K0, _K0P),
    lambda: material_eos.VinetEOS(_RHO0, _K0, _K0P),
    lambda: material_eos.InterpolatedEOS([0.0, 1.0e6, 2.0e6], [5000.0, 4000.0, 3000.0]),
], ids=["constant", "birch_murnaghan", "vinet", "interpolated"])
def test_binary_roundtrip(factory, tmp_path):
    eos = factory()
    rho_before = eos.calc_density(3.0e10, 0.0, 0.5e6)
    reloaded = _binary_copy(eos, make_material_eos(eos.model_name), tmp_path)
    assert reloaded.model_name == eos.model_name
    assert reloaded.calc_density(3.0e10, 0.0, 0.5e6) == pytest.approx(rho_before, rel=1e-12)


def test_isinstance_chain():
    eos = material_eos.BirchMurnaghanEOS(_RHO0, _K0, _K0P)
    assert isinstance(eos, material_eos.MaterialEOSBase)
    assert isinstance(eos, PhysicsBase)
    assert isinstance(eos, TidalPyBaseClass)


# =====================================================================================================================
# Invalid parameters are rejected
# =====================================================================================================================
@pytest.mark.parametrize("model", ("constant", "birch_murnaghan", "vinet"))
@pytest.mark.parametrize("density", (-1000.0, 0.0, math.inf, math.nan))
def test_reference_density_must_be_positive_and_finite(model, density):
    with pytest.raises(ValueError, match="reference density"):
        make_material_eos(model, {"reference_density_kg_m3": density})


@pytest.mark.parametrize("model", ("birch_murnaghan", "vinet"))
@pytest.mark.parametrize("bulk_modulus", (-1.0e11, 0.0, math.inf))
def test_reference_bulk_modulus_must_be_positive_and_finite(model, bulk_modulus):
    with pytest.raises(ValueError, match="reference bulk modulus"):
        make_material_eos(model, {"reference_bulk_modulus_pa": bulk_modulus})


@pytest.mark.parametrize("key", ("shear_modulus_static_pa", "bulk_modulus_static_pa"))
def test_static_moduli_must_not_be_negative(key):
    with pytest.raises(ValueError, match="non-negative"):
        make_material_eos("constant", {key: -1.0})


@pytest.mark.parametrize("call,match", [
    (lambda: make_material_eos("interpolate", {"radius_m": [0.0, 1.0e6], "density_kg_m3": [5000.0, -1.0]}),
     "density"),
    (lambda: material_eos.InterpolatedEOS([0.0, 1.0e6], [5000.0]), None),
    (lambda: make_material_eos("not_a_model"), None),
], ids=["negative_interpolated_density", "interpolated_length_mismatch", "unknown_model"])
def test_invalid_construction_raises(call, match):
    with pytest.raises(ValueError, match=match):
        call()


# =====================================================================================================================
# Thermal terms
# =====================================================================================================================
def test_default_models_are_athermal():
    eos = material_eos.BirchMurnaghanEOS(_RHO0, _K0, _K0P)
    assert eos.thermal_expansion == 0.0
    assert eos.reference_temperature == pytest.approx(_T_REF)
    assert eos.calc_density(1.0e10, 2000.0) == eos.calc_density(1.0e10)


@pytest.mark.parametrize("cls_name", _MODEL_NAMES)
def test_density_at_the_reference_temperature_is_athermal(cls_name):
    """At the reference temperature the density is athermal; hotter is less dense."""
    eos = _thermal_model(cls_name)
    assert eos.thermal_expansion == pytest.approx(_ALPHA)
    athermal = eos.calc_density(2.0e10, None, 0.5e6)
    assert eos.calc_density(2.0e10, _T_REF, 0.5e6) == pytest.approx(athermal, rel=1e-13)
    assert eos.calc_density(2.0e10, 2000.0, 0.5e6) < athermal


@pytest.mark.parametrize("cls_name", ["ConstantDensityEOS", "InterpolatedEOS"])
@pytest.mark.parametrize("temperature", [100.0, 1500.0, 4000.0])
def test_thermal_expansion_factor(cls_name, temperature):
    """Models with no pressure law scale their density by exp(-alpha (T - T_ref))."""
    eos = _thermal_model(cls_name)
    expected = eos.calc_density(0.0, None, 0.5e6) * math.exp(-_ALPHA * (temperature - _T_REF))
    assert eos.calc_density(0.0, temperature, 0.5e6) == pytest.approx(expected, rel=1e-13)


@pytest.mark.parametrize("cls_name,law", _LAWS)
@pytest.mark.parametrize("pressure", [0.0, 1.0e10, 1.0e11])
@pytest.mark.parametrize("temperature", [100.0, 1500.0, 4000.0])
def test_thermal_pressure_roundtrip(cls_name, law, pressure, temperature):
    """The cold law at the solved compression returns the pressure less alpha0 K0 (T - T_ref)."""
    eta = _thermal_model(cls_name).calc_density(pressure, temperature) / _RHO0
    cold_pressure = pressure - _ALPHA * _K0 * (temperature - _T_REF)
    assert getattr(material_eos, law)(eta, _K0, _K0P) == pytest.approx(cold_pressure, rel=1e-6, abs=1.0)


@pytest.mark.parametrize("cls_name", ["BirchMurnaghanEOS", "VinetEOS"])
def test_free_surface_expansion_is_alpha_delta_t(cls_name):
    """At zero pressure a small temperature rise lowers the density by alpha dT to first order."""
    relative_change = 1.0 - _thermal_model(cls_name).calc_density(0.0, _T_REF + 10.0) / _RHO0
    assert relative_change == pytest.approx(_ALPHA * 10.0, rel=1e-2)


@pytest.mark.parametrize("cls_name", ["BirchMurnaghanEOS", "VinetEOS"])
@pytest.mark.parametrize("pressure", [0.0, 1.0e10, 1.0e11])
@pytest.mark.parametrize("temperature", [None, 2500.0])
def test_bulk_modulus_matches_the_density_derivative(cls_name, pressure, temperature):
    """K = rho dP/drho, checked against a centered difference of calc_density."""
    eos = _thermal_model(cls_name)
    step = 1.0e6
    density = eos.calc_density(pressure, temperature)
    slope = eos.calc_density(pressure + step, temperature) - eos.calc_density(pressure - step, temperature)
    assert eos.calc_bulk_modulus(pressure, temperature) == pytest.approx(density * 2.0 * step / slope, rel=1e-6)


@pytest.mark.parametrize("cls_name", ["BirchMurnaghanEOS", "VinetEOS"])
def test_bulk_modulus_at_the_reference_state(cls_name):
    """K is K0 at the reference state, stiffens with compression, and softens with heating."""
    eos = _thermal_model(cls_name)
    assert eos.calc_bulk_modulus(0.0) == pytest.approx(_K0, rel=1e-12)
    assert eos.calc_bulk_modulus(5.0e10) > _K0
    assert eos.calc_bulk_modulus(0.0, 2000.0) < _K0


def test_bulk_modulus_of_models_without_a_pressure_law():
    assert math.isnan(material_eos.ConstantDensityEOS(_RHO0).calc_bulk_modulus(1.0e10, 1000.0))
    table = material_eos.InterpolatedEOS([0.0, 2.0e6], [5000.0, 3000.0], bulk_modulus=[3.0e11, 1.0e11])
    assert table.calc_bulk_modulus(0.0, None, 1.0e6) == pytest.approx(2.0e11)


@pytest.mark.parametrize("cls_name", _MODEL_NAMES)
def test_thermal_terms_survive_config_and_binary_roundtrips(cls_name, tmp_path):
    eos = _thermal_model(cls_name)
    density_before = eos.calc_density(3.0e10, 1800.0, 0.5e6)

    cfg = eos.get_config_dict()
    assert cfg["thermal_expansion_1_k"] == pytest.approx(_ALPHA)
    assert cfg["reference_temperature_k"] == pytest.approx(_T_REF)
    rebuilt = make_material_eos(cfg["model"], {key: value for key, value in cfg.items() if key != "model"})
    assert rebuilt.get_config_dict() == cfg
    assert rebuilt.calc_density(3.0e10, 1800.0, 0.5e6) == pytest.approx(density_before, rel=1e-13)

    loaded = _binary_copy(eos, make_material_eos(eos.model_name), tmp_path)
    assert loaded.thermal_expansion == pytest.approx(_ALPHA)
    assert loaded.calc_density(3.0e10, 1800.0, 0.5e6) == pytest.approx(density_before, rel=1e-13)


def test_reference_temperature_is_configurable():
    eos = make_material_eos("constant", {
        "reference_density_kg_m3": _RHO0,
        "thermal_expansion_1_k": _ALPHA,
        "reference_temperature_k": 1600.0})
    assert eos.reference_temperature == pytest.approx(1600.0)
    assert eos.calc_density(0.0, 1600.0) == pytest.approx(_RHO0, rel=1e-13)
