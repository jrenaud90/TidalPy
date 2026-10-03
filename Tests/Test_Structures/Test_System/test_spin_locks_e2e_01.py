"""Spin-orbit resonance locks end to end: a coupled spin, orbit, and mantle-temperature integration with the locks in a
plain LSODA right-hand side, against the same integration without locks at a tighter tolerance."""
import math

import numpy as np
import pytest

from TidalPy.constants import G, au, year
from TidalPy.Structures import System, build_world

solve_ivp = pytest.importorskip("scipy.integrate").solve_ivp

MYR = 1.0e6 * year
RADIUS = 6.371e6
DENSITY = 5500.0
MASS = (4.0 / 3.0) * math.pi * RADIUS**3 * DENSITY
HEAT_CAPACITY = 1200.0           # [J kg-1 K-1]
CONDUCTIVITY = 3.3               # [W m-1 K-1]
EXPANSIVITY = 3.0e-5             # [K-1]
CRITICAL_RAYLEIGH = 1100.0
SURFACE_TEMPERATURE = 300.0      # [K]
RADIOGENIC_POWER = 2.0e13        # [W]
LOCK_TOLERANCE = 1.0e-5
# A warm Earth-size planet at 0.05 au, e = 0.3, despinning from 10 n: it is captured near 3:2 within the span.
ECCENTRICITY = 0.3
SPIN_RATIO = 10.0
TEMPERATURE = 1700.0             # [K]
SPAN = 2.0                       # [Myr]


def planet_config():
    """A one-layer Andrade mantle with closed-form homogeneous Love numbers and peridotite melting."""
    material = {
        "solid": {
            "thermal_conductivity_w_mk": CONDUCTIVITY,
            "heat_capacity_j_kgk": HEAT_CAPACITY,
            "eos": {"model": "constant", "reference_density_kg_m3": DENSITY, "bulk_modulus_pa": 1.0e12},
            "shear_modulus": {"model": "constant", "shear_modulus_pa": 6.0e10},
            "shear_viscosity": {
                "model": "reference",
                "reference_viscosity_pas": 1.0e21,
                "reference_temperature_k": 1600.0,
                "molar_activation_energy_j_mol": 3.0e5,
                "molar_activation_volume_m3_mol": 0.0,
            },
            "shear_rheology": {"model": "andrade", "alpha": 0.3, "zeta": 1.0},
        },
        "liquid": {"preset": "peridotite"},
        "melting": {"preset": "peridotite"},
        "latent_heat_j_kg": 4.0e5,
    }
    return {
        "name": "planet",
        "type": "terrestrial",
        "radius_m": RADIUS,
        "mass_kg": MASS,
        "eos_solver": {"solve_temperature": False},
        "tides": {"love_method": "homogeneous", "eccentricity_trunc_lvl": 10},
        "layers": {"mantle": {"radius_fraction": 1.0, "temperature_k": 1600.0, "use_melting": True,
                              "material": material}},
    }


class Model:
    """State [a / a0, e, s = spin / n, T / 1000 K], with rates per Myr."""

    def __init__(self):
        self.planet = build_world(planet_config())
        self.system = System("close_in")
        self.system.add_world(build_world("sol"), is_star=True)
        self.system.add_world(self.planet, tidal_host=0, semi_major_axis=0.05 * au, eccentricity=ECCENTRICITY)
        self.a0 = 0.05 * au
        self.calls = 0

    def mantle_cooling(self, temperature):
        """Convective heat loss [W] with Nu = (Ra / Ra_c)^(1/3), at the mantle's own viscosity."""
        delta_t = max(temperature - SURFACE_TEMPERATURE, 0.0)
        viscosity = float(np.asarray(self.planet.get_shear_viscosity(np.array([0.5 * RADIUS])))[0])
        gravity = G * MASS / RADIUS**2
        diffusivity = CONDUCTIVITY / (DENSITY * HEAT_CAPACITY)
        rayleigh = DENSITY * gravity * EXPANSIVITY * delta_t * RADIUS**3 / (diffusivity * viscosity)
        nusselt = max((rayleigh / CRITICAL_RAYLEIGH) ** (1.0 / 3.0), 1.0)
        return 4.0 * math.pi * RADIUS * CONDUCTIVITY * delta_t * nusselt

    def rates(self, state, use_locks):
        self.calls += 1
        a, eccentricity, spin_ratio, temperature = state[0] * self.a0, max(state[1], 0.0), state[2], state[3] * 1.0e3
        self.system.set_semi_major_axis(self.planet, a)
        self.system.set_eccentricity(self.planet, eccentricity)
        self.planet.mantle.temperature = temperature
        self.planet.solve_eos()
        orbital_frequency = self.system.calc_orbital_frequency(self.planet)
        self.planet.set_spin_frequency(spin_ratio * orbital_frequency)
        row = self.system.calc_world_evolution(self.planet, use_locks=use_locks, lock_tolerance=LOCK_TOLERANCE)
        spin_ratio_rate = (row["dspin_dt"] - spin_ratio * row["dn_dt"]) / orbital_frequency
        temperature_rate = (row["tidal_heating"] + RADIOGENIC_POWER - self.mantle_cooling(temperature)) \
            / (MASS * HEAT_CAPACITY)
        return np.array([row["da_dt"] / self.a0, row["de_dt"], spin_ratio_rate, temperature_rate / 1.0e3]) * MYR, row

    def integrate(self, use_locks, rtol):
        result = solve_ivp(
            lambda t, state: self.rates(state, use_locks)[0],
            (0.0, SPAN),
            np.array([1.0, ECCENTRICITY, SPIN_RATIO, TEMPERATURE / 1.0e3]),
            method="LSODA",
            rtol=rtol,
            atol=1.0e-12)
        assert result.success, result.message
        return result


@pytest.fixture(scope="module")
def runs():
    locked_model, free_model = Model(), Model()
    locked = locked_model.integrate(True, 1.0e-7)
    free = free_model.integrate(False, 1.0e-9)
    return locked_model, locked, free_model, free


def test_locks_reproduce_the_free_integration(runs):
    locked_model, locked, _, free = runs
    a_locked, e_locked, s_locked, temperature_locked = locked.y[:, -1]
    a_free, e_free, s_free, temperature_free = free.y[:, -1]
    # The orbit and the interior agree to the integrations' tolerance.
    assert a_locked == pytest.approx(a_free, rel=1.0e-6)
    assert e_locked == pytest.approx(e_free, rel=1.0e-6)
    assert temperature_locked == pytest.approx(temperature_free, rel=1.0e-6)
    # Both end at the warm equilibrium below 3:2; the held spin has relaxed onto it.
    assert abs(s_free - 1.5) < 1.0e-2
    assert abs(s_locked - s_free) < 0.1 * LOCK_TOLERANCE
    # The spin is held at the end of the span, and it moves more slowly than a free spin at the same state.
    rates, row = locked_model.rates(locked.y[:, -1], True)
    free_rates, _ = locked_model.rates(locked.y[:, -1], False)
    assert row["spin_locked"]
    assert abs(rates[2]) < abs(free_rates[2])
