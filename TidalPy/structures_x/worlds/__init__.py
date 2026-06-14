"""TidalPy structures_x.worlds — C++ world class hierarchy.

Exposes the world classes:

- ``BaseWorld``      — identity, orbital/thermal scalars, bulk geometry.
- ``LayeredWorld``   — a world built from an ordered stack of layers.
- ``GasGiantWorld``  — a layered world representing a gas giant.
- ``StarWorld``      — a star (no layers, no EOS); effective temperature/luminosity.
"""

from TidalPy.structures_x.worlds.base import BaseWorld
from TidalPy.structures_x.worlds.layered import LayeredWorld
from TidalPy.structures_x.worlds.gasgiant import GasGiantWorld
from TidalPy.structures_x.worlds.stellar import StarWorld

__all__ = ["BaseWorld", "LayeredWorld", "GasGiantWorld", "StarWorld"]
