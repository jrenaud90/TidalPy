"""TidalPy structures_x.system - the System class linking worlds into an orbiting group.

A :class:`System` holds worlds, each naming its own tidal host (or none) and a two-body orbit about it, and
optionally a star for insolation; it is the container on which orbital evolution is computed.
"""

from TidalPy.structures_x.system.system import System

__all__ = [
    "System",
]
