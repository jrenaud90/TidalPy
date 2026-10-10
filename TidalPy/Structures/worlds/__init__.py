"""TidalPy Structures.worlds: the C++ world class hierarchy.

``BaseWorld`` (an ordered stack of layers, possibly empty, with the identity, orbital and thermal scalars, bulk
geometry, spin, tides, and the EOS and Love solves), ``TerrestrialWorld`` and ``GasGiantWorld`` (the same with their
own world type), and ``StarWorld`` (adds the effective temperature and luminosity). ``world_from_bytes`` rebuilds a
world from its binary record held in memory, as ``BaseWorld.copy`` and unpickling do. ``EOSResult`` is the ``dict``
``BaseWorld.solve_eos`` returns, with a short summary as its ``repr``.
"""

from TidalPy.Structures.worlds.base import BaseWorld, EOSResult, world_from_bytes
from TidalPy.Structures.worlds.terrestrial import TerrestrialWorld
from TidalPy.Structures.worlds.gasgiant import GasGiantWorld
from TidalPy.Structures.worlds.stellar import StarWorld

__all__ = ["BaseWorld", "TerrestrialWorld", "GasGiantWorld", "StarWorld", "EOSResult", "world_from_bytes"]
