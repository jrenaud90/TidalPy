"""Common-task benchmarks that mirror the operations taught in ``Demos/``.

Every task is a zero-argument callable registered with ``@benchmark``. Expensive one-time construction
happens at module scope (or in a ``setup`` callback) so each recorded time reflects the operation under
study rather than its setup. The task set tracks the demos: building worlds and systems, solving the
equation of state, computing tidal heating and Love numbers, and stepping a system's evolution.
"""
from __future__ import annotations

import math
import os
import tempfile

import numpy as np

from harness import benchmark

from TidalPy.constants import G, au, seconds_per_myr
from TidalPy.Structures.configs import build_world, build_system, load_toml
from TidalPy.Structures.worlds import TerrestrialWorld, StarWorld
from TidalPy.Structures.layers import Layer
from TidalPy.Material import Material, Phase
from TidalPy.Viscosity import make_viscosity
from TidalPy.Rheology import Maxwell, Elastic, Andrade
from TidalPy.RadialSolver import radial_solver, homogeneous_love_numbers
from TidalPy.Cooling import make_cooling
from TidalPy.Radiogenics import make_radiogenics
from TidalPy.Tides.classes import make_tide
from TidalPy.Structures.system import System
from TidalPy.Utilities.logging import set_log_level

# A poorly conditioned solve logs a console warning; console I/O inside a timed loop would be charged to the task.
set_log_level("error")

_WORK_DIR = tempfile.mkdtemp(prefix="tidalpy_perf_")


# =====================================================================================================================
# Construction
# =====================================================================================================================
@benchmark("build_world:earth_simple", group="structures", note="2-layer terrestrial from bundled config")
def _build_earth_simple():
    build_world("earth_simple")


@benchmark("build_world:jupiter_simple", group="structures", note="gas giant from bundled config")
def _build_jupiter_simple():
    build_world("jupiter_simple")


@benchmark("build_world:earth_prem", group="structures", note="4-layer PREM-based terrestrial from bundled config")
def _build_earth_prem():
    build_world("earth_prem")


@benchmark("build_world:io", group="structures", note="3-layer Io (core, mantle, asthenosphere) from bundled config")
def _build_io():
    build_world("io")


@benchmark("build_world:sol", group="structures", note="star from bundled config")
def _build_sol():
    build_world("sol")


# The 2-layer Io dict used by the tidal heating task below, as demos P02 and S01 build their worlds.
_IO_DICT = {
    "schema_version": "0.2.0", "name": "Io", "type": "terrestrial",
    "radius_m": 1.8216e6, "mass_kg": 8.9319e22,
    "layers": {"core": {"material": "simple_iron_core", "layer_index": 0,
                        "radius_outer_m": 9.0e5, "use_tides": False},
               "mantle": {"material": "simple_rock", "layer_index": 1,
                          "radius_fraction": 1.0, "use_tides": True}},
}


@benchmark("build_world:from_dict", group="structures", note="2-layer terrestrial from a Python dict")
def _build_from_dict():
    build_world(_IO_DICT)


@benchmark("build_system:sol_system", group="system", note="star + terrestrial + gas giant from bundled config")
def _build_sol_system():
    build_system("sol_system")


# =====================================================================================================================
# Equation of state
# =====================================================================================================================
_prem = build_world("earth_prem")


@benchmark("solve_eos:earth_prem", group="eos", note="integrate the PREM interior")
def _solve_eos_prem():
    _prem.solve_eos()


_earth_simple = build_world("earth_simple")


@benchmark("solve_eos:earth_simple", group="eos", note="integrate the bundled 2-layer Earth interior")
def _solve_eos_earth_simple():
    _earth_simple.solve_eos()


_eos_io = build_world("io")


@benchmark("solve_eos:io", group="eos", note="integrate the bundled 3-layer Io interior")
def _solve_eos_io():
    _eos_io.solve_eos()


@benchmark("build_world+solve_eos:io", group="structures", note="bundled Io built and its interior integrated")
def _build_and_solve_io():
    build_world("io").solve_eos()


# =====================================================================================================================
# Interior profiles
# =====================================================================================================================
_prem.solve_eos()
_PREM_RADII = np.linspace(0.0, _prem.radius, 200)


@benchmark("profiles:earth_prem_200", group="eos", note="density, gravity, and pressure of PREM at 200 radii")
def _profiles_earth_prem():
    _prem.get_density(_PREM_RADII)
    _prem.get_gravity(_PREM_RADII)
    _prem.get_pressure(_PREM_RADII)


# =====================================================================================================================
# World Love numbers
# =====================================================================================================================
_SEMIDIURNAL = 2.0 * math.pi / (12.42 * 3600.0)
_earth_simple.solve_eos()


@benchmark("love_numbers:world_earth_simple", group="radial_solver",
           note="solve_love_numbers on the bundled Earth at the semidiurnal frequency")
def _love_world_earth_simple():
    _earth_simple.solve_love_numbers(frequency=_SEMIDIURNAL)


_eos_io.solve_eos()


@benchmark("love_numbers:world_io", group="radial_solver",
           note="solve_love_numbers on the bundled 3-layer Io at Io's mean motion")
def _love_world_io():
    _eos_io.solve_love_numbers(frequency=4.11e-5)


# =====================================================================================================================
# World configuration, TOML, and binary round trips
# =====================================================================================================================
@benchmark("world_config:get_config_dict", group="structures", note="get_config_dict of the 4-layer PREM Earth")
def _world_get_config_dict():
    _prem.get_config_dict()


_toml_path = os.path.join(_WORK_DIR, "earth_simple.toml")


@benchmark("world_toml:round_trip", group="structures",
           note="save_to_toml + build_world(load_toml) of the bundled 2-layer Earth")
def _world_toml_round_trip():
    _earth_simple.save_to_toml(_toml_path)
    build_world(load_toml(_toml_path))


_world_binary_path = os.path.join(_WORK_DIR, "earth_prem.tpyb")
_world_restored = TerrestrialWorld("placeholder", 1.0, 1.0)


@benchmark("world_binary:round_trip", group="structures",
           note="save_binary + load_binary of the solved 4-layer PREM Earth")
def _world_binary_round_trip():
    _prem.save_binary(_world_binary_path)
    _world_restored.load_binary(_world_binary_path)


# =====================================================================================================================
# Tidal heating (analytic fixed-Q)
# =====================================================================================================================
_M_JUP = 1.898e27
_A_IO = 4.217e8
_N_IO = math.sqrt(G * _M_JUP / _A_IO**3)

_io = build_world({
    "schema_version": "0.2.0", "name": "Io", "type": "terrestrial",
    "radius_m": 1.8216e6, "mass_kg": 8.9319e22, "spin_frequency_rad_s": _N_IO,
    "layers": {"core": {"material": "simple_iron_core", "layer_index": 0,
                        "radius_outer_m": 9.0e5, "use_tides": False},
               "mantle": {"material": "simple_rock", "layer_index": 1,
                          "radius_fraction": 1.0, "use_tides": True}},
})
_io.set_tide_model(make_tide("fixed_q", {"fixed_k": [0.3], "fixed_q": [100.0]}))
_io.set_tide_config(min_degree_l=2, max_degree_l=2, eccentricity_truncation=2, obliquity_truncation=0)
_io.set_spin_frequency(_N_IO)


@benchmark("tidal_heating:fixed_q", group="tides", note="one calc_tides on a fixed-Q world")
def _tidal_heating_fixed_q():
    _io.calc_tides(_N_IO, _N_IO, 0.0041, 0.0, _A_IO, _M_JUP)


# The bundled 3-layer Io with its rheology tide, as demos P09, P10, and P13 use it. Both Love methods at the bundled
# tide settings (degree 2, the default eccentricity truncation).
_bundled_io = build_world("io")
_bundled_io.solve_eos()
_bundled_io_homogeneous = build_world("io")
_bundled_io_homogeneous.solve_eos()
_bundled_io_homogeneous.set_tide_config(love_method="homogeneous")


@benchmark("tidal_heating:io_radial_solver", group="tides",
           note="one calc_tides on the bundled 3-layer Io, radial_solver Love method")
def _tidal_heating_io_radial_solver():
    _bundled_io.calc_tides(_N_IO, _N_IO, 0.0041, 0.0, _A_IO, _M_JUP)
    _bundled_io.get_tidal_heating()


@benchmark("tidal_heating:io_homogeneous", group="tides",
           note="one calc_tides on the bundled 3-layer Io, homogeneous Love method")
def _tidal_heating_io_homogeneous():
    _bundled_io_homogeneous.calc_tides(_N_IO, _N_IO, 0.0041, 0.0, _A_IO, _M_JUP)
    _bundled_io_homogeneous.get_tidal_heating()


# =====================================================================================================================
# Love numbers (radial solver)
# =====================================================================================================================
_R = 1.8216e6
_RHO = 8.9319e22 / (4.0 / 3.0 * math.pi * _R**3)
_SLICES = 50
_mu_scalar = Maxwell().calc_complex_modulus(60.0e9, 1.0e15, _N_IO)


@benchmark("love_numbers:radial_solver", group="radial_solver", note="one homogeneous radial solve for k/h/l")
def _love_radial_solver():
    homogeneous_love_numbers(_R, _RHO, _mu_scalar, _N_IO, num_slices=_SLICES)


@benchmark("love_numbers:prop_matrix", group="radial_solver",
           note="homogeneous static-incompressible propagation-matrix solve for k/h/l")
def _love_prop_matrix():
    homogeneous_love_numbers(_R, _RHO, 60.0e9 + 0.0j, _N_IO, num_slices=_SLICES,
                             layer_is_incompressible=True, love_method='propagation_matrix')


# 3-layer solid / static-liquid / solid planet solved with the shooting method.
_3L_SLICES = 40
_r_core, _r_ocean = 0.35 * _R, 0.55 * _R
_radius_3l = np.concatenate([
    np.linspace(0.0, _r_core, _3L_SLICES),
    np.linspace(_r_core, _r_ocean, _3L_SLICES),
    np.linspace(_r_ocean, _R, _3L_SLICES)])
_density_3l = np.concatenate([
    np.full(_3L_SLICES, 8000.0), np.full(_3L_SLICES, 5000.0), np.full(_3L_SLICES, 3300.0)])
_bulk_3l = np.full(_radius_3l.size, 200.0e9 + 0j)
_shear_3l = np.concatenate([
    np.full(_3L_SLICES, _mu_scalar), np.full(_3L_SLICES, 0.0 + 0.0j), np.full(_3L_SLICES, _mu_scalar)])
_upper_3l = np.array([_r_core, _r_ocean, _R])
_rho_bulk_3l = float(np.sum(_density_3l[1:] * np.diff(_radius_3l**3)) / _R**3)


@benchmark("love_numbers:radial_solver_3layer", group="radial_solver",
           note="solid / static-liquid / solid shooting solve for k/h/l")
def _love_radial_solver_3layer():
    radial_solver(_radius_3l, _density_3l, _bulk_3l, _shear_3l, _N_IO, _rho_bulk_3l,
                  ("solid", "liquid", "solid"), (False, True, False), (False, False, False),
                  _upper_3l, degree_l=2, solve_for=("tidal",))


# =====================================================================================================================
# Rheology
# =====================================================================================================================
_maxwell = Maxwell()


@benchmark("rheology:complex_modulus", group="rheology", number=100000, note="one Maxwell complex-modulus evaluation")
def _rheology_complex_modulus():
    _maxwell.calc_complex_modulus(60.0e9, 1.0e15, _N_IO)


_andrade = Andrade(alpha=0.3, zeta=1.0)
_FREQUENCY_SWEEP = np.logspace(-8, -3, 1000)


@benchmark("rheology:andrade_frequency_sweep_1000", group="rheology",
           note="Andrade complex moduli over 1000 frequencies (vectorized)")
def _rheology_andrade_sweep():
    _andrade.calc_complex_modulus_vectorize_frequency(60.0e9, 1.0e19, _FREQUENCY_SWEEP)


# =====================================================================================================================
# Viscosity, cooling, and radiogenics
# =====================================================================================================================
_reference_viscosity = make_viscosity("reference", {
    "reference_viscosity_pas": 1.0e21, "reference_temperature_k": 1600.0, "molar_activation_energy_j_mol": 3.0e5})


@benchmark("viscosity:reference", group="physics", number=100000, note="one reference-law viscosity evaluation")
def _viscosity_reference():
    _reference_viscosity.calc_viscosity(1500.0, 1.0e9)


# Demo P10's silicate mantle: delta T, thickness, gravity, density, viscosity, conductivity, diffusivity, expansivity.
_convection = make_cooling("convection")
_VISCOSITY_SWEEP = np.logspace(16.0, 26.0, 200)
_COOLING_ARGS = (1000.0, 9.1e5, 1.5, 3300.0, _VISCOSITY_SWEEP, 3.75, 3.75 / (3300.0 * 1200.0), 5.2e-5)


@benchmark("cooling:convection_viscosity_sweep_200", group="physics",
           note="convective cooling over 200 viscosities (vectorized)")
def _cooling_convection_sweep():
    _convection.calc_cooling_vectorize_viscosity(*_COOLING_ARGS)


_isotopes = make_radiogenics("isotope", {"isotopes": "modern_day_chondritic"})
_RADIOGENIC_TIMES = np.linspace(0.0, 4600.0 * seconds_per_myr, 10000)


@benchmark("radiogenics:isotope_time_sweep_10k", group="physics",
           note="chondritic isotope heating over 10k times (vectorized)")
def _radiogenics_isotope_sweep():
    _isotopes.calc_heating_vectorize_time(_RADIOGENIC_TIMES, 1.0)


# =====================================================================================================================
# System evolution
# =====================================================================================================================
_star_host = build_world("jupiter_simple")
_evo_system = System("evo")
_evo_system.add_world(_star_host)
_evo_io = build_world({
    "schema_version": "0.2.0", "name": "EvoIo", "type": "terrestrial",
    "radius_m": 1.8216e6, "mass_kg": 8.9319e22, "spin_frequency_rad_s": _N_IO,
    "layers": {"mantle": {"material": "simple_rock", "layer_index": 0,
                          "radius_fraction": 1.0, "use_tides": True}},
})
_evo_io.set_tide_model(make_tide("fixed_q", {"fixed_k": [0.3], "fixed_q": [100.0]}))
_evo_io.set_tide_config(min_degree_l=2, max_degree_l=2, eccentricity_truncation=2, obliquity_truncation=0)
_evo_system.add_world(_evo_io, tidal_host=_star_host, semi_major_axis=_A_IO, eccentricity=0.01)
_evo_io.set_spin_frequency(_N_IO)


@benchmark("system_evolution:fixed_q", group="system", note="one calc_world_evolution step")
def _system_evolution():
    _evo_system.calc_world_evolution(_evo_io)


# The Jupiter-Io system of the migration guide: the bundled 3-layer Io with its rheology tide and radial solves.
_jovian = System("jovian")
_jovian_host = build_world("jupiter_simple")
_jovian_io = build_world("io")
_jovian_io.solve_eos()
_jovian.add_world(_jovian_host)
_jovian.add_world(_jovian_io, tidal_host=_jovian_host, semi_major_axis=_A_IO, eccentricity=0.0041)
_jovian_io.set_spin_frequency(_jovian.calc_orbital_frequency(_jovian_io))


@benchmark("system_evolution:io_rheology", group="system",
           note="calc_world_evolution of the bundled 3-layer Io about Jupiter (rheology tide)")
def _system_evolution_io_rheology():
    _jovian.calc_world_evolution(_jovian_io)


# Demo S03: the bundled Sun, Earth, and Moon, with the Earth and Moon raising tides on each other.
_ems = System("Earth-Moon-Sun")
_ems_earth = build_world("earth_simple")
_ems_moon = build_world("luna")
_ems_earth.solve_eos()
_ems_moon.solve_eos()
_ems.add_world(build_world("sol"), is_star=True)
_ems.add_world(_ems_earth)
_ems.add_world(_ems_moon, tidal_host=_ems_earth, semi_major_axis=3.84748e8, eccentricity=0.0549)
_ems.set_tidal_host(_ems_earth, _ems_moon)
_ems_moon.set_spin_frequency(_ems.calc_orbital_frequency(_ems_moon))
for _world in (_ems_earth, _ems_moon):
    _ems.set_stellar_semi_major_axis(_world, au)
    _ems.set_stellar_eccentricity(_world, 0.0167)


@benchmark("system_evolution:earth_moon_pair", group="system",
           note="calc_pair_evolution of the bundled Earth and Moon (both rheology tides)")
def _system_evolution_earth_moon_pair():
    _ems.calc_pair_evolution(_ems_moon)


@benchmark("system:insolation", group="system", note="insolation flux and equilibrium temperature of the Moon")
def _system_insolation():
    _ems.calc_insolation_flux(_ems_moon)
    _ems.calc_equilibrium_temperature(_ems_moon)


# Demo S01: a star with a terrestrial planet and a gas giant, both fixed-Q, evolved together.
def _fixed_q_planet(config):
    world = build_world(config)
    world.set_tide_model(make_tide("fixed_q", {"fixed_k": [0.3], "fixed_q": [100.0]}))
    world.set_tide_config(min_degree_l=2, max_degree_l=2, eccentricity_truncation=2, obliquity_truncation=0)
    return world


_multi_star = StarWorld("Sun", 6.957e8, 1.988e30)
_multi_star.set_effective_temperature(5772.0)
_multi_inner = _fixed_q_planet({**_IO_DICT, "name": "Inner", "radius_m": 6.0e6, "mass_kg": 5.0e24,
                                "spin_frequency_rad_s": 2.0e-5})
_multi_outer = _fixed_q_planet({
    "schema_version": "0.2.0", "name": "Outer", "type": "gasgiant",
    "radius_m": 6.0e7, "mass_kg": 6.0e26, "spin_frequency_rad_s": 1.0e-4,
    "layers": {"envelope": {"material": "simple_gas", "layer_index": 0,
                            "radius_fraction": 1.0, "use_tides": True}}})
_multi_system = System("multi")
_multi_system.add_world(_multi_star, is_star=True)
_multi_system.add_world(_multi_inner, tidal_host=_multi_star, semi_major_axis=0.20 * au, eccentricity=0.05)
_multi_system.add_world(_multi_outer, tidal_host=_multi_star, semi_major_axis=1.50 * au, eccentricity=0.02)


@benchmark("system_evolution:multi_world", group="system",
           note="calc_system_evolution of a star, a terrestrial planet, and a gas giant (fixed-Q)")
def _system_evolution_multi_world():
    _multi_system.calc_system_evolution()


# =====================================================================================================================
# Equation of state (Birch-Murnaghan)
# =====================================================================================================================
_bm_world = TerrestrialWorld("bm_planet", 6.371e6, 5.972e24)
_bm_material = Material(solid=Phase(
    eos={"model": "birch_murnaghan", "reference_density_kg_m3": 4000.0},
    shear_modulus={"model": "constant", "shear_modulus_pa": 80.0e9},
    shear_viscosity={"model": "constant", "reference_viscosity_pas": 1.0e21},
    bulk_viscosity={"model": "constant", "reference_viscosity_pas": 1.0e30}))
_bm_world.add_layer(Layer(
    "mantle",
    0,
    0.0,
    6.371e6,
    5.972e24,
    _bm_material,
    shear_rheology=Maxwell(),
    bulk_rheology=Elastic()))


@benchmark("solve_eos:birch_murnaghan", group="eos", note="homogeneous terrestrial with a Birch-Murnaghan EOS")
def _solve_eos_birch_murnaghan():
    _bm_world.solve_eos()


# =====================================================================================================================
# 3D tides (rheology tide, fully collapsed total)
# =====================================================================================================================
# A uniform Maxwell mantle (60 GPa, 1e15 Pa s) shared by the 3D tide tasks.
_MAXWELL_MANTLE = Material(solid=Phase(
    eos={"model": "constant", "reference_density_kg_m3": _RHO, "bulk_modulus_pa": 200.0e9},
    shear_modulus={"model": "constant", "shear_modulus_pa": 60.0e9},
    shear_viscosity={"model": "constant", "reference_viscosity_pas": 1.0e15},
    bulk_viscosity={"model": "constant", "reference_viscosity_pas": 1.0e15}))


def _maxwell_mantle_layer():
    """The dynamic Maxwell mantle filling a world of Io's radius and mass."""
    return Layer(
        "mantle",
        0,
        0.0,
        _R,
        8.9319e22,
        _MAXWELL_MANTLE,
        is_static=False,
        shear_rheology=Maxwell(),
        bulk_rheology=Elastic())


_rheo_io = TerrestrialWorld("RheoIo", _R, 8.9319e22)
_rheo_io.add_layer(_maxwell_mantle_layer())
_rheo_io.set_tide_model(make_tide("rheology"))
_rheo_io.set_tide_config(min_degree_l=2, max_degree_l=2, eccentricity_truncation=2, obliquity_truncation=0)
_rheo_io.solve_eos()


@benchmark("tides_3d:collapse_total", group="tides", note="fully collapsed 3D heating total (rheology tide)")
def _tides_3d_collapse_total():
    _rheo_io.calc_3d_tides(_N_IO, 1.5 * _N_IO, 0.0041, 0.0, _A_IO, _M_JUP,
                           latitude_summed=True, longitude_summed=True, radial_summed=True)


# Demo P06: secular heating along one line of points (120 colatitudes just below the surface).
_HEATING_COLATITUDES = np.linspace(0.02, np.pi - 0.02, 120)
_HEATING_RADII = np.full(_HEATING_COLATITUDES.size, 0.99 * _R)


@benchmark("tides_3d:heating_array_120", group="tides",
           note="get_3d_tidal_heating_array at 120 surface colatitudes (rheology tide)")
def _tides_3d_heating_array():
    _rheo_io.get_3d_tidal_heating_array(_N_IO, _N_IO, 0.0041, 0.0, _A_IO, _M_JUP,
                                        _HEATING_RADII, _HEATING_COLATITUDES)


# =====================================================================================================================
# 3D grids on one thread and on every logical core
# =====================================================================================================================
# Degrees 2 to 3 with eccentricity, a non-synchronous spin, and obliquity, so hundreds of waves reach every point. The
# radial solves are the same in both variants; only the per-point evaluation after them runs on the extra threads.
_grid_io = TerrestrialWorld("GridIo", _R, 8.9319e22)
_grid_io.add_layer(_maxwell_mantle_layer())
_grid_io.set_tide_model(make_tide("rheology"))
_grid_io.set_tide_config(min_degree_l=2, max_degree_l=3, eccentricity_truncation=10, obliquity_truncation=2)
_grid_io.solve_eos()
_GRID_STATE = (_N_IO, 1.2 * _N_IO, 0.1, 0.2, _A_IO, _M_JUP)
_GRID_AXES = dict(radii=np.linspace(0.05 * _R, 0.99 * _R, 20),
                  colatitudes=np.linspace(0.02, np.pi - 0.02, 45),
                  longitudes=np.linspace(0.0, 2.0 * np.pi, 90, endpoint=False))
_GRID_TIMES = np.linspace(0.0, 2.0 * np.pi / _N_IO, 24, endpoint=False)
_ALL_THREADS = os.cpu_count() or 1


def _secular_map(num_threads):
    _grid_io.calc_3d_tides(*_GRID_STATE, num_threads=num_threads, **_GRID_AXES)


def _stress_strain(num_threads):
    _grid_io.calc_3d_stress_strain(*_GRID_STATE, times=_GRID_TIMES[:4], num_threads=num_threads, **_GRID_AXES)


def _displacements(num_threads):
    _grid_io.calc_3d_displacements(*_GRID_STATE, times=_GRID_TIMES, num_threads=num_threads, **_GRID_AXES)


@benchmark("tides_3d:secular_map_1_thread", group="tides", repeats=3,
           note="secular heating map 20 x 45 x 90, degrees 2 to 3, e^10, obliquity, 1 thread")
def _tides_3d_secular_map_1_thread():
    _secular_map(1)


@benchmark("tides_3d:secular_map_all_threads", group="tides", repeats=3,
           note=f"secular heating map 20 x 45 x 90, degrees 2 to 3, e^10, obliquity, {_ALL_THREADS} threads")
def _tides_3d_secular_map_all_threads():
    _secular_map(_ALL_THREADS)


@benchmark("tides_3d:stress_strain_1_thread", group="tides", repeats=3,
           note="stress and strain 20 x 45 x 90 x 4 times, degrees 2 to 3, e^10, obliquity, 1 thread")
def _tides_3d_stress_strain_1_thread():
    _stress_strain(1)


@benchmark("tides_3d:stress_strain_all_threads", group="tides", repeats=3,
           note=f"stress and strain 20 x 45 x 90 x 4 times, degrees 2 to 3, e^10, obliquity, {_ALL_THREADS} threads")
def _tides_3d_stress_strain_all_threads():
    _stress_strain(_ALL_THREADS)


@benchmark("tides_3d:displacements_1_thread", group="tides", repeats=3,
           note="displacements 20 x 45 x 90 x 24 times, degrees 2 to 3, e^10, obliquity, 1 thread")
def _tides_3d_displacements_1_thread():
    _displacements(1)


@benchmark("tides_3d:displacements_all_threads", group="tides", repeats=3,
           note=f"displacements 20 x 45 x 90 x 24 times, degrees 2 to 3, e^10, obliquity, {_ALL_THREADS} threads")
def _tides_3d_displacements_all_threads():
    _displacements(_ALL_THREADS)


# =====================================================================================================================
# System binary round trip
# =====================================================================================================================
_sol_system = build_system("sol_system")
_binary_path = os.path.join(_WORK_DIR, "sol_system_roundtrip.tpb")


@benchmark("system_binary:round_trip", group="system", note="save_binary + load_binary of the bundled sol system")
def _system_binary_round_trip():
    _sol_system.save_binary(_binary_path)
    System("loaded").load_binary(_binary_path)
