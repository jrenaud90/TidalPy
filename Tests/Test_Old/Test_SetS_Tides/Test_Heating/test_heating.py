import TidalPy


import numpy as np

from TidalPy.tides.heating import calculate_volumetric_heating

def test_volumetric_heating():
    """ Test the volumetric heating function. """

    # Ensure imports work
    from TidalPy.tides.heating import calculate_volumetric_heating
    from TidalPy.tides import calculate_volumetric_heating

    stress = (10. + 0.1j) * np.ones((6, 10, 12), dtype=np.complex128)
    strain = (5. - 0.1j) * np.ones((6, 10, 12), dtype=np.complex128)

    volumetric_heating = calculate_volumetric_heating(stress, strain)

    assert volumetric_heating.shape == (10, 12)
    assert volumetric_heating.dtype == np.float64
    assert np.all(volumetric_heating >= 0.)
