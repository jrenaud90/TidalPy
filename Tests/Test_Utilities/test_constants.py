"""Tests for loading and updating the constants in ``TidalPy.constants``."""
import math

import pytest
import scipy.constants


@pytest.mark.parametrize(
    "name, expected",
    [
        ("radius_solar", 6.957e8),
        ("mass_solar", 1.988435e30),
        ("mass_earth", 5.9721986e24),
        ("radius_earth", 6.371008e6),
        ("ppm", 1.e-6),
        ("ppb", 1.e-9),
        # Config-loaded values: these hold only for an unmodified default config.
        ("min_frequency", 1.0e-16),
        ("max_frequency", 1.0e8),
        ("G", scipy.constants.G),
        ("R", scipy.constants.R),
        ("k_boltzman", scipy.constants.k),
        ("sbc", scipy.constants.Stefan_Boltzmann),
    ])
def test_loading_constants(name, expected):
    """Compile-time, config-loaded, and SciPy-sourced constants load with their expected values."""
    import TidalPy.constants

    assert math.isclose(getattr(TidalPy.constants, name), expected)


def test_update_constants():
    """A config update changes a static constant, and a second update restores it."""
    import TidalPy

    from TidalPy.constants import test_constant
    assert math.isclose(test_constant, 42.0)
    del test_constant

    TidalPy.reinit(dict(numerical=dict(test_constant=9001.0)))

    from TidalPy.constants import test_constant
    assert math.isclose(test_constant, 9001.0)
    del test_constant

    TidalPy.reinit(dict(numerical=dict(test_constant=42.0)))

    from TidalPy.constants import test_constant
    assert math.isclose(test_constant, 42.0)
