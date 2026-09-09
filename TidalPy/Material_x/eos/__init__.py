"""TidalPy.Material_x.eos — C++ equation-of-state models and solver.

Exposes the material EOS model hierarchy (density as a function of pressure or
radius) and its name-based factory:

- ``ConstantDensityEOS`` (alias ``"constant"``) — incompressible.
- ``BirchMurnaghanEOS``  (alias ``"bm"``)       — 3rd-order Birch-Murnaghan.
- ``VinetEOS``           (alias ``"vinet"``)    — Vinet/UBER EOS.
- ``InterpolatedEOS``    (alias ``"interp"``)   — density(radius) lookup table.
"""

from TidalPy.Material_x.eos.material_eos import (
    MaterialEOSBase,
    ConstantDensityEOS,
    BirchMurnaghanEOS,
    VinetEOS,
    InterpolatedEOS,
    make_material_eos,
    birch_murnaghan_pressure,
    vinet_pressure,
)

__all__ = [
    "MaterialEOSBase",
    "ConstantDensityEOS",
    "BirchMurnaghanEOS",
    "VinetEOS",
    "InterpolatedEOS",
    "make_material_eos",
    "birch_murnaghan_pressure",
    "vinet_pressure",
]
