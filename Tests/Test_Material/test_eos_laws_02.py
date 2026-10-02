"""Equation-of-state laws at the edges of their range: the density inversion in deep tension and compression, where a
pressure law turns over, and tabulated profiles read on uneven tables, before and after a binary round trip."""

import math

import numpy as np
import pytest

from TidalPy.Material.laws import make_eos, make_shear_modulus
from TidalPy.Viscosity import make_viscosity

_RHO0 = 3500.0
_K0 = 1.30e11


def _bm_pressure(eta, k0, kp):
    return 1.5 * k0 * (eta ** (7.0 / 3.0) - eta ** (5.0 / 3.0)) * (1.0 + 0.75 * (kp - 4.0) * (eta ** (2.0 / 3.0) - 1.0))


def _vinet_pressure(eta, k0, kp):
    x = eta ** (-1.0 / 3.0)
    return 3.0 * k0 * (1.0 - x) / x ** 2 * math.exp(1.5 * (kp - 1.0) * (1.0 - x))


_LAWS = {"birch_murnaghan": _bm_pressure, "vinet": _vinet_pressure}


def _law(name, k0_prime, k0=_K0):
    return make_eos(name, {"reference_density_kg_m3": _RHO0, "reference_bulk_modulus_pa": k0,
                           "bulk_modulus_derivative": k0_prime})


def _binary_copy(model, tmp_path):
    path = str(tmp_path / "law.tpyb")
    model.save_binary(path)
    blank = type(model)()
    blank.load_binary(path)
    return blank


@pytest.mark.parametrize("name", sorted(_LAWS))
@pytest.mark.parametrize("k0_prime", [3.2, 4.0, 5.5])
@pytest.mark.parametrize("pressure", [-5.0e9, -1.0e6, 1.0e3, 1.0e9, 1.4e11, 2.5e11])
def test_inversion_is_tight_in_compression_and_tension(name, k0_prime, pressure):
    """Inside the monotonic range the inversion recovers the pressure to its tolerance."""
    eta = _law(name, k0_prime).calc_density(pressure) / _RHO0
    assert (eta > 1.0) == (pressure > 0.0)
    # K is the slope, so a compression good to 1e-13 gives a pressure good to about K 1e-13.
    assert _LAWS[name](eta, _K0, k0_prime) == pytest.approx(pressure, rel=1.0e-10, abs=_K0 * 1.0e-11)


@pytest.mark.parametrize("name", sorted(_LAWS))
@pytest.mark.parametrize("k0_prime", [3.2, 4.0, 5.5])
def test_density_is_continuous_where_the_law_turns_over_in_tension(name, k0_prime):
    """Past the law's minimum pressure the density holds the turning point, with no jump: the structure solve probes
    deep tension while its central pressure is a guess, and a jump would stall the stepper."""
    law = _law(name, k0_prime)
    floor_density = law.calc_density(-10.0 * _K0)
    assert law.calc_density(-100.0 * _K0) == floor_density
    floor_eta = floor_density / _RHO0
    assert 0.0 < floor_eta < 1.0
    floor_pressure = _LAWS[name](floor_eta, _K0, k0_prime)
    assert _LAWS[name](0.99 * floor_eta, _K0, k0_prime) > floor_pressure
    # The law is flat at its minimum, so 1e-9 K0 in pressure is about its square root in compression.
    just_inside = law.calc_density(floor_pressure + 1.0e-9 * _K0)
    assert just_inside == pytest.approx(floor_density, rel=1.0e-3)
    assert just_inside >= floor_density


def test_birch_murnaghan_holds_its_maximum_pressure_when_it_turns_over_in_compression():
    """With K0' < 4 the law turns over at large compression and the density holds; with K0' >= 4 it keeps rising."""
    law = _law("birch_murnaghan", 3.2)
    ceiling_density = law.calc_density(1.0e15)
    assert law.calc_density(1.0e16) == ceiling_density
    ceiling_eta = ceiling_density / _RHO0
    assert _bm_pressure(1.01 * ceiling_eta, _K0, 3.2) < _bm_pressure(ceiling_eta, _K0, 3.2)
    unbounded = _law("birch_murnaghan", 4.5)
    assert unbounded.calc_density(1.0e16) > unbounded.calc_density(1.0e15)


@pytest.mark.parametrize("name", sorted(_LAWS))
def test_a_reloaded_law_inverts_like_the_one_that_was_saved(name, tmp_path):
    """The monotonic range is derived from K0 and K0', so a binary load must rebuild it."""
    law = _law(name, 3.4, k0=2.2e11)
    reloaded = _binary_copy(law, tmp_path)
    for pressure in (-1.0e13, -1.0e9, 5.0e10, 1.0e15):
        assert reloaded.calc_density(pressure) == law.calc_density(pressure)


def test_birch_murnaghan_and_vinet_agree_at_small_compression():
    assert _law("birch_murnaghan", 4.5).calc_density(5.0e9) == pytest.approx(
        _law("vinet", 4.5).calc_density(5.0e9), rel=0.02)


@pytest.mark.parametrize("law", ["polytrope", "constant"])
def test_a_law_without_a_compression_scaling_keeps_alpha0(law):
    """The polytrope (a barotrope with no reference density) and the constant law (whose density scales by
    exp(-alpha0 (T - T_ref))) report alpha0 at any pressure."""
    assert make_eos(law, {"thermal_expansion_1_k": 3.0e-5}).calc_eos(1.0e11, 2000.0)["thermal_expansion"] == 3.0e-5


def test_tabulated_laws_read_as_numpy_on_an_uneven_table(tmp_path):
    """A PREM-like table (uneven, with a repeated radius at a discontinuity) reads as numpy.interp in every tabulated
    law, before and after a binary round trip."""
    radius = np.sort(np.concatenate([
        np.linspace(0.0, 3.0e6, 7), [3.0e6], np.linspace(3.2e6, 6.0e6, 9),
        6.0e6 + np.cumsum(np.geomspace(1.0e5, 2.0e3, 40))]))
    rows = np.arange(radius.size, dtype=float)
    tables = {"density": 1.0e4 - 10.0 * rows, "bulk_modulus": 1.0e11 + 3.0e8 * rows,
              "shear_modulus": 5.0e10 + 1.0e8 * rows ** 1.5, "viscosity": 10.0 ** (18.0 + 0.05 * rows)}
    # Distinct values across the discontinuity, so reading the wrong side of it would show.
    discontinuity = int(np.flatnonzero(np.diff(radius) == 0.0)[0])
    for values in tables.values():
        values[discontinuity + 1:] *= 0.7
    eos = make_eos("interpolate", {"radius_m": radius, "density_kg_m3": tables["density"],
                                   "bulk_modulus_pa": tables["bulk_modulus"]})
    shear = make_shear_modulus("interpolate", {"radius_m": radius, "shear_modulus_pa": tables["shear_modulus"]})
    viscosity = make_viscosity("interpolate", {"radius_m": radius, "viscosity_pas": tables["viscosity"]})
    queries = np.concatenate([radius, radius[:-1] + 0.5 * np.diff(radius), np.nextafter(radius, -np.inf),
                              np.linspace(-1.0e6, radius[-1] + 1.0e6, 3001)])

    def check(eos_law, shear_law, viscosity_law):
        state = eos_law.calc_eos(0.0, radius=queries)
        reads = {"density": state["density"], "bulk_modulus": state["bulk_modulus"],
                 "shear_modulus": shear_law.calc_shear_modulus(0.0, radius=queries),
                 "viscosity": viscosity_law.calc_viscosity(1000.0, 0.0, queries)}
        for name, got in reads.items():
            np.testing.assert_allclose(got, np.interp(queries, radius, tables[name]), rtol=1.0e-14, atol=0.0,
                                       err_msg=name)

    check(eos, shear, viscosity)
    check(_binary_copy(eos, tmp_path), _binary_copy(shear, tmp_path), _binary_copy(viscosity, tmp_path))
