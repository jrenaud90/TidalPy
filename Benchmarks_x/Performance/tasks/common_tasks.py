"""Common-task benchmarks that mirror the operations taught in ``Demos_x/``.

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

from TidalPy.constants import G
from TidalPy.structures_x.configs import build_world, build_system
from TidalPy.structures_x.worlds.layered import LayeredWorld
from TidalPy.structures_x.layers.physics import PhysicsLayer
from TidalPy.Material_x.eos.material_eos import ConstantDensityEOS, BirchMurnaghanEOS
from TidalPy.viscosity_x import make_viscosity
from TidalPy.rheology_x import Maxwell, Elastic
from TidalPy.RadialSolver_x import radial_solver, homogeneous_love_numbers
from TidalPy.Tides_x.classes import make_tide
from TidalPy.structures_x.system import System


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


# =====================================================================================================================
# Tidal heating (analytic fixed-Q)
# =====================================================================================================================
_M_JUP = 1.898e27
_A_IO = 4.217e8
_N_IO = math.sqrt(G * _M_JUP / _A_IO**3)

_io = build_world({
    "schema_version": "0.2.0", "name": "Io", "type": "terrestrial",
    "radius_m": 1.8216e6, "mass_kg": 8.9319e22, "spin_frequency_rad_s": _N_IO,
    "layers": {"core": {"class": "physics", "type": "iron", "layer_index": 0,
                        "radius_outer_m": 9.0e5, "is_tidal": False},
               "mantle": {"class": "solidliquid", "type": "mantle_rock", "layer_index": 1,
                          "radius_fraction": 1.0, "is_tidal": True}},
})
_io.set_tide_model(make_tide("fixed_q", {"fixed_k": [0.3], "fixed_q": [100.0]}))
_io.set_tide_config(min_degree_l=2, max_degree_l=2, eccentricity_truncation=2, obliquity_truncation=0)
_io.set_spin_frequency(_N_IO)


@benchmark("tidal_heating:fixed_q", group="tides", note="one calc_tides on a fixed-Q world")
def _tidal_heating_fixed_q():
    _io.calc_tides(_N_IO, _N_IO, 0.0041, 0.0, _A_IO, _M_JUP)


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
                             layer_is_incompressible=True, use_prop_matrix=True)


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


# =====================================================================================================================
# System evolution
# =====================================================================================================================
_star_host = build_world("jupiter_simple")
_evo_system = System("evo")
_evo_system.add_world(_star_host, is_host=True)
_evo_io = build_world({
    "schema_version": "0.2.0", "name": "EvoIo", "type": "terrestrial",
    "radius_m": 1.8216e6, "mass_kg": 8.9319e22, "spin_frequency_rad_s": _N_IO,
    "layers": {"mantle": {"class": "solidliquid", "type": "mantle_rock", "layer_index": 0,
                          "radius_fraction": 1.0, "is_tidal": True}},
})
_evo_io.set_tide_model(make_tide("fixed_q", {"fixed_k": [0.3], "fixed_q": [100.0]}))
_evo_io.set_tide_config(min_degree_l=2, max_degree_l=2, eccentricity_truncation=2, obliquity_truncation=0)
_evo_system.add_world(_evo_io, semi_major_axis=_A_IO, eccentricity=0.01)
_evo_io.set_spin_frequency(_N_IO)


@benchmark("system_evolution:fixed_q", group="system", note="one calc_world_evolution step")
def _system_evolution():
    _evo_system.calc_world_evolution(_evo_io)


# =====================================================================================================================
# Equation of state (Birch-Murnaghan)
# =====================================================================================================================
_bm_world = LayeredWorld("bm_planet", 6.371e6, 5.972e24)
_bm_layer = PhysicsLayer("mantle", 0, 0.0, 6.371e6, 5.972e24,
                         shear_modulus_static_pa=80.0e9, bulk_modulus_static_pa=200.0e9)
_bm_layer.set_eos(BirchMurnaghanEOS(reference_density_kg_m3=4000.0))
_bm_layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity": 1.0e21}))
_bm_layer.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity": 1.0e30}))
_bm_layer.set_shear_rheology(Maxwell())
_bm_layer.set_bulk_rheology(Elastic())
_bm_world.add_layer(_bm_layer)


@benchmark("solve_eos:birch_murnaghan", group="eos", note="homogeneous terrestrial with a Birch-Murnaghan EOS")
def _solve_eos_birch_murnaghan():
    _bm_world.solve_eos()


# =====================================================================================================================
# 3D tides (rheology tide, fully collapsed total)
# =====================================================================================================================
_rheo_io = LayeredWorld("RheoIo", _R, 8.9319e22)
_rheo_layer = PhysicsLayer("mantle", 0, 0.0, _R, 8.9319e22,
                           shear_modulus_static_pa=60.0e9, bulk_modulus_static_pa=200.0e9)
_rheo_layer.is_static = False
_rheo_layer.set_eos(ConstantDensityEOS(reference_density_kg_m3=_RHO))
_rheo_layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity": 1.0e15}))
_rheo_layer.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity": 1.0e15}))
_rheo_layer.set_shear_rheology(Maxwell())
_rheo_layer.set_bulk_rheology(Elastic())
_rheo_io.add_layer(_rheo_layer)
_rheo_io.set_tide_model(make_tide("rheology"))
_rheo_io.set_tide_config(min_degree_l=2, max_degree_l=2, eccentricity_truncation=2, obliquity_truncation=0)
_rheo_io.solve_eos()


@benchmark("tides_3d:collapse_total", group="tides", note="fully collapsed 3D heating total (rheology tide)")
def _tides_3d_collapse_total():
    _rheo_io.calc_3d_tides(_N_IO, 1.5 * _N_IO, 0.0041, 0.0, _A_IO, _M_JUP,
                           latitude_summed=True, longitude_summed=True, radial_summed=True)


# =====================================================================================================================
# System binary round trip
# =====================================================================================================================
_sol_system = build_system("sol_system")
_binary_path = os.path.join(tempfile.mkdtemp(prefix="tidalpy_perf_"), "sol_system_roundtrip.tpb")


@benchmark("system_binary:round_trip", group="system", note="save_binary + load_binary of the bundled sol system")
def _system_binary_round_trip():
    _sol_system.save_binary(_binary_path)
    System("loaded").load_binary(_binary_path)
