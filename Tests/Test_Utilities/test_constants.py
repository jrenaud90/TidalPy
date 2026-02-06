import pytest
import math

def test_loading_constants():
    """Test that constants are able to be loaded from TidalPy and do a spot check for some values."""

    # True constants (defined by static constexpr in "constants_.hpp").
    from TidalPy.constants import radius_solar, mass_solar, radius_earth, mass_earth, ppm, ppb

    assert math.isclose(radius_solar, 6.957e8)
    assert math.isclose(mass_solar, 1.988435e30)
    assert math.isclose(mass_earth, 5.9721986e24)
    assert math.isclose(radius_earth, 6.371008e6)
    assert math.isclose(ppm, 1.e-6)
    assert math.isclose(ppb, 1.e-9)

    # Test static members that are loaded from configuration file. Note that this test will only work on TidalPy that
    #  has an unmodified config from default.
    from TidalPy.constants import min_frequency, max_frequency, min_spin_orbit_diff

    assert math.isclose(min_frequency, 1.0e-14)
    assert math.isclose(max_frequency, 1.0e8)
    assert math.isclose(min_spin_orbit_diff, 1.0e-10)

    # Test constants loaded from 3rd party packages.
    from scipy.constants import G as sp_G, R as sp_R, k as sp_k, Stefan_Boltzmann as sp_sbc
    from TidalPy.constants import G, R, k_boltzman, sbc

    assert math.isclose(sp_G, G)
    assert math.isclose(sp_R, R)
    assert math.isclose(sp_k, k_boltzman)
    assert math.isclose(sp_sbc, sbc)


def test_update_constants():
    """Test the ability to change static constants and for that information to not be lost."""

    import TidalPy

    from TidalPy.constants import test_constant
    assert math.isclose(test_constant, 42.0)
    del test_constant
    
    new_config = dict(
        debug = dict(
            test_constant = 9001.0
        )
    )
    TidalPy.reinit(new_config)

    from TidalPy.constants import test_constant
    assert math.isclose(test_constant, 9001.0)
    del test_constant

    new_config = dict(
        debug = dict(
            test_constant = 42.0
        )
    )
    TidalPy.reinit(new_config)

    from TidalPy.constants import test_constant
    assert math.isclose(test_constant, 42.0)
