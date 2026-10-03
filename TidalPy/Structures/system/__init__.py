"""TidalPy Structures.system - the System class linking worlds into an orbiting group.

A :class:`System` holds worlds, each naming its own tidal host (or none) and a two-body orbit about it, and
optionally a star for insolation; it is the container on which orbital evolution is computed. ``system_from_bytes``
rebuilds a system from its binary record held in memory, as ``System.copy`` and unpickling do.
"""

from TidalPy.Structures.system.system import System, system_from_bytes

__all__ = [
    "System",
    "system_from_bytes",
]
