"""TidalPy structures_x.system - the System class linking worlds into an orbiting group.

A :class:`System` holds a host world plus orbiting worlds, each with a two-body orbit about the host,
and is the container on which orbital evolution is computed.
"""

from TidalPy.structures_x.system.system import System

__all__ = [
    "System",
]
