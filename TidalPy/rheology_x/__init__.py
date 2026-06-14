"""TidalPy rheology_x — C++ rheology model hierarchy.

Exposes the seven rheology models and a name-based factory:

- ``Elastic``  (alias ``"off"``)
- ``Viscous``  (alias ``"newton"``)
- ``Voigt``    (alias ``"voigt-kelvin"``)
- ``Maxwell``
- ``Burgers``
- ``Andrade``
- ``Sundberg`` (alias ``"sundberg-cooper"``)

Each model computes the complex modulus (shear μ* or bulk K*) [Pa].
"""

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
