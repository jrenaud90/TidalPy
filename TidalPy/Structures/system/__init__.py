"""TidalPy Structures.system - the System class linking worlds into an orbiting group.

A :class:`System` holds worlds, each naming its own tidal host (or none) and a two-body orbit about it, and
optionally a star for insolation; it is the container on which orbital evolution is computed, and
:meth:`System.evolve` integrates a pair of its worlds (:class:`PairedEvolutionResult`, one
:class:`EvolutionResult` per world). ``system_from_bytes`` rebuilds a system from its
binary record held in memory, as ``System.copy`` and unpickling do.
"""

from TidalPy.Structures.system.system import (
    EvolutionResult,
    PairedEvolutionResult,
    System,
    build_spin_window,
    evolution_result_from_state,
    locate_spin_root,
    paired_evolution_result_from_state,
    refine_spin_root,
    spin_search_points,
    system_from_bytes,
)

__all__ = [
    "EvolutionResult",
    "PairedEvolutionResult",
    "System",
    "build_spin_window",
    "evolution_result_from_state",
    "locate_spin_root",
    "paired_evolution_result_from_state",
    "refine_spin_root",
    "spin_search_points",
    "system_from_bytes",
]
