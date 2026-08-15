"""Common-task benchmarks that mirror the operations taught in ``Demos_x/``.

Every task is a zero-argument callable registered with ``@benchmark``. Keep the timed call itself
minimal: do expensive one-time construction in a module-level fixture or a ``setup`` callback so the
recorded time reflects the operation under study rather than its setup.

The task set grows alongside the demos so the two stay in step. As of now it covers world and system
construction from the bundled configs; tidal, Love-number, and evolution tasks are added as those
demos land.
"""
from __future__ import annotations

from harness import benchmark

from TidalPy.structures_x.configs import build_world, build_system


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
