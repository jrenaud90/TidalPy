def test_exoplanet_download():
    """Tests the exoplanet archive download utility."""

    from TidalPy.utilities.exoplanets import get_exoplanet_data

    exoplanet_data = get_exoplanet_data(
        ensure_radius=True,
        radius_cutoff=1.5,
        ensure_mass=True,
        mass_cutoff=None,
        ensure_orbital_period=True,
        orbital_period_cutoff=50,
        ensure_eccentricity=True,
        eccentricity_threshold=None,
        star_type=None,
        only_defaults=True)
    
    assert len(exoplanet_data) > 0

def test_exoplanet_download_stellar_type():
    """Tests the exoplanet archive download utility using the different kinds of stellar types."""

    from TidalPy.utilities.exoplanets import get_exoplanet_data

    # Try just a single stellar type
    exoplanet_data = get_exoplanet_data(
        ensure_radius=True,
        radius_cutoff=1.5,
        ensure_mass=True,
        mass_cutoff=None,
        ensure_orbital_period=True,
        orbital_period_cutoff=50,
        ensure_eccentricity=True,
        eccentricity_threshold=None,
        star_type='M',
        only_defaults=True)
    
    assert len(exoplanet_data) > 0

    # Try just a list of stellar types
    exoplanet_data = get_exoplanet_data(
        ensure_radius=True,
        radius_cutoff=1.5,
        ensure_mass=True,
        mass_cutoff=None,
        ensure_orbital_period=True,
        orbital_period_cutoff=50,
        ensure_eccentricity=True,
        eccentricity_threshold=None,
        star_type=['F', 'G', 'K'],
        only_defaults=True)
    
    assert len(exoplanet_data) > 0
