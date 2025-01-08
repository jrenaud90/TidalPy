from time import time

import astropy
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive

from TidalPy.logger import get_logger

log = get_logger("TidalPy")


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

    # Build criteria string
    criteria = []

    if ensure_radius:
        criteria.append('pl_rade>0')

    if radius_cutoff is not None:
        criteria.append(f'pl_rade<={radius_cutoff}')
    
    if ensure_mass:
        criteria.append('(pl_masse>0 OR pl_msinie>0)')
    
    if mass_cutoff is not None:
        criteria.append(f'(pl_masse<={mass_cutoff} OR pl_msinie<={mass_cutoff})')
    
    if ensure_orbital_period:
        criteria.append('pl_orbper>0')
    
    if orbital_period_cutoff is not None:
        criteria.append(f'pl_orbper<={orbital_period_cutoff}')

    if ensure_eccentricity:
        criteria.append(f'pl_orbeccen>0')
    
    if eccentricity_threshold:
        criteria.append(f'pl_orbeccen>={eccentricity_threshold}')

    if star_type is not None:
        # The exoplanet archive uses % as a wildcard
        if type(star_type) in (list, tuple):
            # Create a list of stellar types
            star_criteria = ''
            num_types = len(star_type)
            for stype_i, stype in enumerate(star_type):
                if stype_i < (num_types - 1):
                    star_criteria += f"st_spectype like '{stype}%' OR "
                else:
                    # Last item, no "or"
                    star_criteria += f"st_spectype like '{stype}%'"
        else:
            star_criteria = f"st_spectype like '{star_type}%'"
        
        criteria.append(f'({star_criteria})')
    
    if only_defaults:
        criteria.append('default_flag=1')
    
    if len(criteria) > 0:
        criteria_str = " AND ".join(criteria)
    else:
        criteria_str = ""

    log.debug(f'Getting exoplanet data using the criteria: {criteria}.')

    t0 = time()
    try:
        exoplanet_data = NasaExoplanetArchive.query_criteria(
            table="ps",
            select="*",
            where=criteria_str
            )
    except Exception as QueryError:
        log.exception(QueryError)
        raise QueryError
    log.debug(f'Exoplanet data retrieved in {time() - t0: 0.2f} seconds.')
    num_exoplanets = len(exoplanet_data)
    log.debug(f'{num_exoplanets} exoplanets were found which met the criteria.')

    return exoplanet_data

