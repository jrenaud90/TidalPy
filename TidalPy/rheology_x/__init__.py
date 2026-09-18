"""TidalPy rheology_x: the C++ rheology model hierarchy.

Each model returns the complex modulus (shear μ* or bulk K*) [Pa]: ``Elastic`` (alias ``"off"``),
``Viscous`` (alias ``"newton"``), ``Voigt`` (alias ``"voigt-kelvin"``), ``Maxwell``, ``Burgers``,
``Andrade``, and ``Sundberg`` (alias ``"sundberg-cooper"``), plus the ``make_rheology`` factory.
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
