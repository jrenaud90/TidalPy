"""C++ rheology models and their name-based factory. Each returns a complex modulus (shear or bulk) [Pa]."""

from TidalPy.rheology_x.rheology import (
    RheologyBase,
    Elastic,
    Viscous,
    Voigt,
    Maxwell,
    Burgers,
    Andrade,
    Sundberg,
    make_rheology,
    elastic,
    viscous,
    voigt,
    maxwell,
    burgers,
    andrade,
    sundberg,
)

__all__ = [
    # Model classes
    "RheologyBase",
    "Elastic",
    "Viscous",
    "Voigt",
    "Maxwell",
    "Burgers",
    "Andrade",
    "Sundberg",
    # Factory
    "make_rheology",
    # Direct complex-modulus convenience functions (float or np.ndarray inputs)
    "elastic",
    "viscous",
    "voigt",
    "maxwell",
    "burgers",
    "andrade",
    "sundberg",
]
