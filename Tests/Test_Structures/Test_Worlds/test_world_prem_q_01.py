"""PREM with its own attenuation: a radial profile's quality factors, not a viscosity, set the loss (q_provided)."""

import copy
import math

import numpy as np
import pytest

from TidalPy.constants import G
from TidalPy.Rheology import seismic_q
from TidalPy.Structures import build_world
from TidalPy.Structures.configs import data_file, worldpack

_M2 = 2.0 * math.pi / (12.4206012 * 3600.0)   # the principal lunar semidiurnal tide [rad/s]
_PREM_REFERENCE = 2.0 * math.pi                # PREM's 1 s reference period [rad/s]


@pytest.fixture(autouse=True)
def _packaged_worldpack(tmp_path, monkeypatch):
    """Resolve bundled files from the package, not from a data-directory copy an older install may have left."""
    monkeypatch.setattr(worldpack, "get_worlds_dir", lambda: str(tmp_path))


def _prem_q_config(**changes):
    config = copy.deepcopy(build_world("earth_prem_q").portable_config)
    config.update(changes)
    return config


def _k2(world, frequency=_M2):
    world.solve_eos(G_to_use=G, verbose=False)
    result = world.solve_love_numbers(frequency=frequency, degree_l=2)
    assert result["success"] is True, result["message"]
    return complex(world.love_number_k)


def _prem_mapping(with_q=True, with_viscosity=False):
    """PREM as a mapping of arrays, the way a profile reaches build_world from Python."""
    arrays = data_file.load_radial_data(worldpack.resolve_data_file("PREM.csv"))
    mapping = {"radius_m": arrays["radius_m"], "density_kg_m3": arrays["density_kg_m3"],
               "vp_m_s": arrays["vp_m_s"], "vs_m_s": arrays["vs_m_s"]}
    if with_q:
        mapping.update(q_mu=arrays["shear_q"], q_kappa=arrays["bulk_q"])
    if with_viscosity:
        mapping["shear_viscosity_pas"] = np.full_like(arrays["radius_m"], 1.0e21)
    return mapping


def _mapping_world(name, mapping):
    """A PREM-like world config whose profile is a mapping of arrays, with its quality factors switched on."""
    return {"schema_version": "0.2.0", "name": name, "type": "terrestrial", "radius_m": 6.371e6,
            "mass_kg": 5.972e24, "data": mapping, "q_provided": True}


# =====================================================================================================================
# The profile's quality factors
# =====================================================================================================================
def test_prem_file_carries_prems_quality_factors():
    arrays = data_file.load_radial_data(worldpack.resolve_data_file("PREM.csv"))
    radius_km = arrays["radius_m"] / 1.0e3

    def q_at(radius):
        return float(arrays["shear_q"][np.argmin(np.abs(radius_km - radius))])

    # Dziewonski and Anderson (1981): crust and lid 600, low-velocity zone 80, upper mantle 143, lower mantle 312,
    # outer core 0 (liquid), inner core 84.6; Q_kappa 57823 except 1327.7 in the inner core.
    assert [q_at(r) for r in (6300.0, 6200.0, 6000.0, 4000.0, 2500.0, 500.0)] == [600.0, 80.0, 143.0, 312.0, 0.0,
                                                                                    84.6]
    assert set(np.unique(arrays["bulk_q"])) == {1327.7, 57823.0}
    assert np.all((arrays["shear_q"] == 0.0) == (arrays["shear_modulus_pa"] == 0.0))


def test_quality_factor_columns_are_ignored_without_q_provided():
    """earth_prem reads the same file and stays elastic."""
    config = build_world("earth_prem").source_config
    for layer in config["layers"].values():
        assert "shear_viscosity_pas" not in layer["material"]
        assert "shear_rheology" not in layer
    assert abs(_k2(build_world("earth_prem")).imag) < 1.0e-12


# =====================================================================================================================
# The built world
# =====================================================================================================================
def test_earth_prem_q_gives_each_solid_layer_seismic_q():
    world = build_world("earth_prem_q")
    live = world.get_config_dict()["layers"]
    solid = [name for name, layer in zip(live, world) if layer.is_solid]
    assert solid == ["layer_0", "layer_2"]
    for name in solid:
        for table in ("shear_rheology", "bulk_rheology"):
            assert live[name][table] == {"model": "seismic_q", "reference_frequency_rad_s": _PREM_REFERENCE,
                                         "q_frequency_exponent": 0.0}
    mantle_q = world.source_config["layers"]["layer_2"]["material"]["shear_viscosity_pas"]
    assert (min(mantle_q), max(mantle_q)) == (80.0, 600.0)
    # The liquid outer core has no quality factor and no rheology.
    outer_core = world.source_config["layers"]["layer_1"]
    assert "shear_viscosity_pas" not in outer_core["material"]
    assert "shear_rheology" not in outer_core


def test_mantle_modulus_is_the_seismic_q_of_the_profile():
    """The complex modulus the Love solve sees is seismic_q applied to the profile's modulus and Q."""
    world = build_world("earth_prem_q")
    world.solve_eos(G_to_use=G, verbose=False)
    mantle = list(world)[2]
    radius = 5.0e6                                    # lower mantle, Q_mu = 312
    shear = mantle.get_shear_modulus(radius)
    got = mantle.calc_complex_shear_modulus(radius, _M2)
    expected = seismic_q(shear, 312.0, _M2)
    assert got.real == pytest.approx(expected.real, rel=1.0e-12)
    assert got.imag == pytest.approx(expected.imag, rel=1.0e-12)
    assert got.real / got.imag == pytest.approx(312.0, rel=1.0e-12)


def test_at_the_reference_frequency_k2_is_the_elastic_prems():
    """PREM's moduli are its 1 s moduli, so at 1 s only the loss is new."""
    elastic = _k2(build_world("earth_prem"), _PREM_REFERENCE)
    lossy = _k2(build_world("earth_prem_q"), _PREM_REFERENCE)
    assert lossy.real == pytest.approx(elastic.real, rel=1.0e-4)
    assert lossy.imag < 0.0


def test_dispersion_softens_the_tidal_response():
    """At M2 the dispersion raises Re k2 above the elastic value, and a Q falling with period adds more loss."""
    elastic = _k2(build_world("earth_prem"))
    constant_q = _k2(build_world("earth_prem_q"))
    power_law = _k2(build_world(_prem_q_config(q_frequency_exponent=0.15)))
    assert elastic.real < constant_q.real < power_law.real
    assert 0.0 > constant_q.imag > power_law.imag
    # The effective Q of k2 at M2: about 500 with PREM's Q alone, about 100 with a = 0.15.
    assert abs(constant_q) / abs(constant_q.imag) == pytest.approx(495.0, rel=0.02)
    assert abs(power_law) / abs(power_law.imag) == pytest.approx(100.0, rel=0.02)


def test_reference_frequency_moves_the_dispersion():
    """The same Q quoted at a lower reference frequency sits closer to M2, so M2 is stiffer."""
    at_1s = _k2(build_world("earth_prem_q"))
    at_200s = _k2(build_world(_prem_q_config(q_reference_frequency_rad_s=2.0 * math.pi / 200.0)))
    assert at_200s.real < at_1s.real


def test_a_layer_table_overrides_the_world_settings():
    config = _prem_q_config()
    config["layers"] = {"layer_0": {"shear_rheology": {"model": "constant_q", "q_frequency_exponent": 0.2}},
                        "mantle": {"layer_index": 2, "bulk_rheology": {"model": "elastic"}}}
    live = build_world(config).get_config_dict()["layers"]
    assert live["layer_0"]["shear_rheology"]["q_frequency_exponent"] == 0.2
    assert live["layer_0"]["shear_rheology"]["reference_frequency_rad_s"] == _PREM_REFERENCE
    assert live["layer_0"]["bulk_rheology"]["model"] == "seismic_q"
    assert live["mantle"]["bulk_rheology"] == {"model": "elastic"}
    assert live["mantle"]["shear_rheology"]["model"] == "seismic_q"


def test_a_mapping_profile_takes_q_provided_too():
    config = _mapping_world("PREM-Q-mapping", _prem_mapping())
    assert _k2(build_world(config)) == pytest.approx(_k2(build_world("earth_prem_q")), rel=1.0e-12)


def test_saved_world_keeps_q_provided(tmp_path):
    world = build_world("earth_prem_q")
    path = str(tmp_path / "prem_q.toml")
    world.save_to_toml(path)
    reloaded = build_world(path)
    assert reloaded.portable_config["q_provided"] is True
    assert _k2(reloaded) == pytest.approx(_k2(world), rel=1.0e-12)


# =====================================================================================================================
# Refusals
# =====================================================================================================================
@pytest.mark.parametrize("changes,match", [
    ({"q_provided": "yes"}, "must be true or false"),
    ({"q_provided": False}, "not 'q_provided = true'"),
    ({"q_frequency_exponent": 1.0}, "q_frequency_exponent"),
    ({"q_frequency_exponent": -0.1}, "q_frequency_exponent"),
    ({"q_reference_frequency_rad_s": 0.0}, "q_reference_frequency_rad_s"),
])
def test_bad_q_settings_are_refused(changes, match):
    with pytest.raises(ValueError, match=match):
        build_world(_prem_q_config(**changes))


def test_q_settings_need_a_profile():
    config = copy.deepcopy(build_world("earth_simple").get_config_dict())
    config["q_provided"] = True
    with pytest.raises(ValueError, match="gives no profile"):
        build_world(config)


@pytest.mark.parametrize("mapping,match", [
    (dict(with_q=False), "no Q_mu column"),
    (dict(with_viscosity=True), "gives viscosities too"),
])
def test_profile_must_give_q_and_no_viscosity(mapping, match):
    config = _mapping_world("PREM-Q-bad", _prem_mapping(**mapping))
    with pytest.raises(ValueError, match=match):
        build_world(config)


@pytest.mark.parametrize("layer_table,match", [
    ({"shear_rheology": {"model": "maxwell"}}, "only the 'seismic_q' rheology reads"),
    ({"bulk_rheology": {"model": "andrade"}}, "or 'elastic', which ignores them"),
    ({"material": {"shear_viscosity_static_pas": 1.0e21}}, "not a viscosity"),
    ({"type": "rock"}, "Leave 'type' unset"),
    # A solid layer's viscosity slot holds Q_mu, which a melt law or a Rayleigh number would read as a viscosity.
    ({"material": {"partial_melt": {"model": "henning"}}}, "Only 'off' is allowed"),
    ({"material": {"partial_melt": {"model": "spohn"}}}, "Only 'off' is allowed"),
    ({"cooling": {"model": "convective"}}, "Use 'conduction' or 'off'"),
])
def test_layer_tables_that_would_misread_q_are_refused(layer_table, match):
    config = _prem_q_config(layers={"layer_2": layer_table})
    with pytest.raises(ValueError, match=match):
        build_world(config)


def test_models_that_do_not_read_the_viscosity_are_allowed():
    layer_table = {"material": {"partial_melt": {"model": "off"}}, "cooling": {"model": "conduction"}}
    world = build_world(_prem_q_config(layers={"layer_2": layer_table}))
    assert world.layer_2.get_config_dict()["cooling"]["model"] == "conduction"


def test_a_solid_layer_with_zero_q_is_refused():
    mapping = _prem_mapping()
    mapping["q_mu"] = np.where(mapping["radius_m"] > 6.0e6, 0.0, mapping["q_mu"])
    config = _mapping_world("PREM-Q-zero", mapping)
    with pytest.raises(ValueError, match="quality factors must be positive"):
        build_world(config)
