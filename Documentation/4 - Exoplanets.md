# TidalPy and Exoplanet Data
TidalPy has built in functionality to download data from NASA's exoplanet archive (using the `astroquery` package). You must be connected to the internet and may run into issues if you have a firewall or strict permissions on your computer / network.

## Downloading Exoplanet Data

Below is the function description with details on its arguments.

```python
def get_exoplanet_data(
        ensure_radius: bool = True,
        radius_cutoff: float = 1.5,
        ensure_mass: bool = True,
        mass_cutoff: float = None,
        ensure_orbital_period: bool  = True,
        orbital_period_cutoff: float = 50,
        ensure_eccentricity: float = True,
        eccentricity_threshold: float = None,
        star_type: str = None,
        only_defaults: bool = True) -> astropy.table.table.QTable:
    """Utility to download exoplanet data from NASA's exoplanet archive.

    Utilizes the `astroquery` package and NEXSCI's exoplanet archive data and API.

    Parameters
    ----------
    ensure_radius : bool, optional
        Require the planet(s) to have a measured radius, by default True.
    radius_cutoff : float, optional
        Require the planet(s) to have a radius <= to this value, by default 1.5. [R_Earth]
        If set to None then no cutoff will be used.
    ensure_mass : bool, optional
        Require the planet(s) to have a measured mass or minimum mass, by default True.
    mass_cutoff : float, optional
        Require the planet(s) to have a mass or minimum mass <= to this value, by default None. [M_Earth]
        If set to None then no cutoff will be used.
    ensure_orbital_period : bool, optional
        Require the planet(s) to have a measured orbital period, by default True.
    orbital_period_cutoff : float, optional
        Require the planet(s) to have an orbital period <= to this value, by default 50. [days]
    ensure_eccentricity : float, optional
        Require the planet(s) to have a measured non-zero eccentricity, by default True.
    eccentricity_threshold : float, optional
        Require the planet(s) to have an eccentricity >= to this value, by default None.
    star_type : str, optional
        Restrict to planets orbiting this stellar type, by default None.
        For example, `star_type='M'` or `star_type=['K', 'G', 'K']`
    only_defaults : bool, optional
        Only return the default values as set by the NASA Exoplanet Archive, by default True.
        If this is set to False then you may get multiple copies of the same planet with different parameters.

    Returns
    -------
    exoplanet_data : astropy.table.table.QTable
        An astropy QTable with rows for each exoplanet found that meets the criteria and columns for all available 
        data. Where applicable, data will have dimensional units utilizing `astropy.units.core.Unit` classes.

    Raises
    ------
    QueryError
        An error occurred during the NASA Exoplanet Archive query, this could have been a network connection problem,
        a permissions problem, a server issue, or that the criteria used was badly formatted.
    """
```

Below is how to run the function.

```python
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

print(f'{len(exoplanet_data)} exoplanets found matching the criteria.')

print(exoplanet_data)
```
