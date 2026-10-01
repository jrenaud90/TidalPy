"""The Anderson-Gruneisen thermal expansivity of the material EOS models."""
import math

import pytest

from TidalPy.Material.eos.material_eos import BirchMurnaghanEOS, ConstantDensityEOS, InterpolatedEOS, make_material_eos


def test_constant_expansivity_is_the_default():
    """Without an Anderson-Gruneisen parameter the expansivity is alpha0 at every density."""
    model = BirchMurnaghanEOS(reference_density=3300.0, thermal_expansion=3.0e-5)
    assert model.anderson_gruneisen_parameter == 0.0
    for density in (2000.0, 3300.0, 6000.0):
        assert model.calc_thermal_expansion(density) == 3.0e-5


@pytest.mark.parametrize("delta, kappa", [(5.5, 0.0), (5.5, 1.4), (3.0, 2.0)])
def test_expansivity_follows_the_anderson_gruneisen_closed_form(delta, kappa):
    """alpha = alpha0 exp[(d0 / k) ((rho0 / rho)^k - 1)], and alpha0 (rho0 / rho)^d0 for k = 0."""
    model = BirchMurnaghanEOS(reference_density=3300.0, thermal_expansion=3.0e-5,
                              anderson_gruneisen_parameter=delta, anderson_gruneisen_exponent=kappa)
    for density in (3300.0, 4500.0, 5500.0):
        ratio = 3300.0 / density
        expected = 3.0e-5 * (ratio ** delta if kappa == 0.0 else math.exp((delta / kappa) * (ratio ** kappa - 1.0)))
        assert model.calc_thermal_expansion(density) == pytest.approx(expected, rel=1e-14)
    assert model.calc_thermal_expansion(3300.0) == pytest.approx(3.0e-5, rel=1e-15)


def test_a_falling_delta_softens_the_drop():
    """A delta_T that falls with compression (kappa > 0) leaves the compressed expansivity higher than a constant one."""
    constant = BirchMurnaghanEOS(reference_density=3300.0, thermal_expansion=3.0e-5, anderson_gruneisen_parameter=5.5)
    falling = BirchMurnaghanEOS(reference_density=3300.0, thermal_expansion=3.0e-5, anderson_gruneisen_parameter=5.5,
                                anderson_gruneisen_exponent=1.4)
    assert falling.calc_thermal_expansion(5500.0) > constant.calc_thermal_expansion(5500.0)


def test_models_without_a_reference_density_keep_alpha0():
    """The interpolated profile has no reference density, so its expansivity stays alpha0."""
    model = InterpolatedEOS([0.0, 1.0e6], [5000.0, 4000.0], thermal_expansion=2.0e-5,
                            anderson_gruneisen_parameter=5.0)
    assert model.calc_thermal_expansion(4500.0) == 2.0e-5


def test_parameters_round_trip_through_the_config():
    """get_config_dict writes both parameters and the factory reads them back."""
    model = ConstantDensityEOS(reference_density=3300.0, thermal_expansion=3.0e-5, anderson_gruneisen_parameter=5.5,
                               anderson_gruneisen_exponent=1.4)
    config = model.get_config_dict()
    assert config["anderson_gruneisen_parameter"] == 5.5
    assert config["anderson_gruneisen_exponent"] == 1.4
    rebuilt = make_material_eos(config.pop("model"), config)
    assert rebuilt.anderson_gruneisen_parameter == 5.5
    assert rebuilt.anderson_gruneisen_exponent == 1.4
    assert rebuilt.calc_thermal_expansion(4000.0) == model.calc_thermal_expansion(4000.0)
