"""Common-task benchmarks that mirror the operations taught in ``Demos_x/``.

Every task is a zero-argument callable registered with ``@benchmark``. Expensive one-time construction
happens at module scope (or in a ``setup`` callback) so each recorded time reflects the operation under
study rather than its setup. The task set tracks the demos: building worlds and systems, solving the
equation of state, computing tidal heating and Love numbers, and stepping a system's evolution.
"""
from __future__ import annotations

import math

import numpy as np

from harness import benchmark

from TidalPy.constants import G
from TidalPy.structures_x.configs import build_world, build_system
from TidalPy.structures_x.worlds.layered import LayeredWorld
from TidalPy.structures_x.layers.physics import PhysicsLayer
from TidalPy.Material_x.eos.material_eos import ConstantDensityEOS
from TidalPy.viscosity_x import make_viscosity
from TidalPy.rheology_x import Maxwell
from TidalPy.RadialSolver_x import radial_solver
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
_radius = np.linspace(0.0, _R, _SLICES)
_density = np.full(_SLICES, _RHO)
_complex_bulk = np.full(_SLICES, 200.0e9 + 0j)
_complex_shear = np.full(_SLICES, Maxwell().calc_complex_modulus(60.0e9, 1.0e15, _N_IO), dtype=np.complex128)
_upper = np.array([_R])


@benchmark("love_numbers:radial_solver", group="radial_solver", note="one homogeneous radial solve for k/h/l")
def _love_radial_solver():
    radial_solver(_radius, _density, _complex_bulk, _complex_shear, _N_IO, _RHO,
                  ("solid",), (True,), (False,), _upper, degree_l=2, solve_for=("tidal",))


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
