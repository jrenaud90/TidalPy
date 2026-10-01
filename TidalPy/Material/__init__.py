"""Materials: the property laws (equations of state, shear-modulus laws), and the phases and materials that combine
them with viscosity, rheology, and melting laws into a material's state at a point."""

from TidalPy.Material.laws import (
    make_eos,
    make_shear_modulus,
)
from TidalPy.Material.material import (
    Phase,
    Material,
    make_phase,
    make_material,
)

__all__ = [
    "make_eos",
    "make_shear_modulus",
    "Phase",
    "Material",
    "make_phase",
    "make_material",
]
