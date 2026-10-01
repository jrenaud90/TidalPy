"""Materials: the property laws (equations of state, shear-modulus laws), the phases and materials that combine
them with viscosity, rheology, and melting laws into a material's state at a point, and MatPack, the named materials
TidalPy ships."""

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
from TidalPy.Material.matpack import (
    CATEGORIES,
    available_materials,
    install_matpack,
    load_material,
    material_config,
    material_info,
    merge_material_tables,
)

__all__ = [
    "make_eos",
    "make_shear_modulus",
    "Phase",
    "Material",
    "make_phase",
    "make_material",
    "CATEGORIES",
    "available_materials",
    "install_matpack",
    "load_material",
    "material_config",
    "material_info",
    "merge_material_tables",
]
