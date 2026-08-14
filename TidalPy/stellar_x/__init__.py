"""TidalPy stellar_x - C++ stellar physics model hierarchy.

Exposes the three luminosity models and a name-based factory:

- ``FixedLuminosity``    (alias ``"constant"``)              - luminosity set directly.
- ``MassToLuminosity``   (aliases ``"cuntz_wang"``/``"cw"``) - piecewise main-sequence L(M).
- ``PowerLawLuminosity`` (alias ``"power_law"``)             - single power law L ~ M^p.

Each model computes a star's luminosity [W] from its mass via ``calc_luminosity`` and shares the
Stefan-Boltzmann effective-temperature conversions.
"""

from TidalPy.stellar_x.luminosity import (
    LuminosityBase,
    FixedLuminosity,
    MassToLuminosity,
    PowerLawLuminosity,
    make_luminosity,
    fixed,
    mass_to_luminosity,
    power_law,
)

__all__ = [
    # Model classes
    "LuminosityBase",
    "FixedLuminosity",
    "MassToLuminosity",
    "PowerLawLuminosity",
    # Factory
    "make_luminosity",
    # Direct luminosity convenience functions (float or np.ndarray mass)
    "fixed",
    "mass_to_luminosity",
    "power_law",
]
