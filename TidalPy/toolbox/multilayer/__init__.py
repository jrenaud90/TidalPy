from .homogen_solid import calculate_homogen_solid
from .liquid_solid import calculate_ls
from .solid_liquid_solid import calculate_sls
from .solid_solid_liquid_solid import calculate_ssls

KNOWN_INTERIOR_MODELS = {
    'homogeneous': calculate_homogen_solid,
    'liquid-solid': calculate_ls,
    'solid-liquid-solid': calculate_sls,
    'solid-solid-liquid-solid': calculate_ssls
    }
