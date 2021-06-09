"""quick_tides.py module.

Purpose
-------
The primary purpose of this module is to quickly and easily provide the user with tidal dissipation calculation
    functions (heating and orbit-spin derivatives) for homogeneous worlds without requiring the user to manually
    pull together several different TidalPy functions.

Limitations
-----------
These functions may have an impact to performance and certainly to flexibility. For example, there are various
    checks performed in each that a user may want to skip. It is recommended to only use this function as a starting
    point and then develop a custom function or script to suite your specific needs.

These are intended to be used for homogeneous worlds with one layer. Dissipation for multi-layer worlds can be done
    with TidalPy's OOP approach or a custom script.

"""

import copy
from typing import Tuple, Dict, TYPE_CHECKING, Union

import numpy as np

from .conversions import days2rads, orbital_motion2semi_a
from ..dynamics import semia_eccen_derivatives_dual, semia_eccen_derivatives, spin_rate_derivative
from ..exceptions import MissingArgumentError, ArgumentException, IncorrectArgumentType, MissingAttributeError
from ..rheology.complex_compliance import known_models as known_compliance_models
from ..rheology.complex_compliance.complex_compliance import compliance_dict_helper
from ..tides.ctl_funcs import linear_dt
from ..tides.dissipation import calc_tidal_susceptibility
from ..tides.methods.global_approx import ctl_neg_imk_helper_func, cpl_neg_imk_helper_func
from ..tides.mode_manipulation import find_mode_manipulators
from ..utilities.types import FloatArray, ComplexArray, NoneType

if TYPE_CHECKING:
    from structures.world_types import all_tidal_world_types

NoneFloatArray = Union[NoneType, FloatArray]
NoneFloat = Union[NoneType, float]
NoneInt = Union[NoneType, int]
SingleBodyResultType = Dict[str, Union[FloatArray, Dict[int, Union[ComplexArray, FloatArray]]]]


def quick_tidal_dissipation(host_mass: float, target_radius: float, target_mass: float,
                            target_gravity: float, target_density: float, target_moi: float,
                            viscosity: FloatArray = None, shear_modulus: FloatArray = None, rheology: str = 'Maxwell',
                            complex_compliance_inputs: Tuple[float, ...] = None,
                            eccentricity: FloatArray = None, obliquity: FloatArray = None,
                            orbital_frequency: FloatArray = None, orbital_period: FloatArray = None,
                            spin_frequency: FloatArray = None, spin_period: FloatArray = None,
                            max_tidal_order_l: Union[int, NoneType] = 2,
                            eccentricity_truncation_lvl: Union[int, NoneType] = 2,
                            use_obliquity: Union[bool, NoneType] = True,
                            tidal_scale: float = 1., fixed_k2: float = 0.3, fixed_q: float = 100.,
                            fixed_dt: float = None,
                            dspin_dt_scale: float = 1., de_dt_scale: float = 1., da_dt_scale: float = 1.,
                            precalculated_mode_results: dict = None,
                            calculate_orbit_spin_derivatives: bool = False) -> SingleBodyResultType:
    """ Calculate the single body tidal dissipation for a target world orbiting its tidal host.

    Assumes the host's dissipation is zero (no contribution to the orbital evolution derivatives).

    This function pulls together several TidalPy packages to offer an easy to use interface to calculate tidal
        dissipation. However it may come with an impact to performance and flexibility. For example, there are various
        checks performed that may want to skip. It is recommended to only use this function as a starting point and
        then develop a custom function or script to suite your specific needs.

    Parameters
    ----------
    host_mass : float
    target_radius : float
    target_mass : float
    target_gravity : float
    target_density : float
    target_moi : float
    viscosity : FloatArray
    shear_modulus : FloatArray
    rheology : str
    complex_compliance_inputs : Tuple[float, ...]
    eccentricity : FloatArray
    obliquity : FloatArray
    orbital_frequency : FloatArray
    orbital_period : FloatArray
    spin_frequency : FloatArray
    spin_period : FloatArray
    max_tidal_order_l : int = 2
    eccentricity_truncation_lvl : int = 2,
    use_obliquity : bool = True
    tidal_scale : float = 1.
    fixed_k2 : float = 0.3
    fixed_q : float = 100.
    fixed_dt : float = None
        If `None`, a fixed_dt will be calculated as 1 / (Q * n)
    dspin_dt_scale : float = 1.
    de_dt_scale : float = 1.
    da_dt_scale : float = 1.
    precalculated_mode_results : dict = None
        Several things can be pre-calculated to save on performance if doing, for example, dual body analysis
    calculate_orbit_spin_derivatives: bool = False
        If `True`, then the function will calculate the spin and orbital derivatives assuming the
            host is non-dissipative.

    Returns
    -------
    dissipation_results : dict
        Dissipation results for the target world
        'tidal_heating': FloatArray,
        'dUdM': FloatArray,
        'dUdw': FloatArray,
        'dUdO': FloatArray,
        'love_number_by_orderl': List[ComplexArray],
        'negative_imk_by_orderl': List[FloatArray],
        'effective_q_by_orderl': List[FloatArray]
    """

    # Flag that tracks if array version of TidalPy's functions should be used or not.
    use_array = False

    # Get orbital frequency
    if orbital_frequency is None:
        if orbital_period is None:
            raise MissingArgumentError('Orbital period or frequency is required.')
        else:
            orbital_frequency = days2rads(orbital_period)
    semi_major_axis = orbital_motion2semi_a(orbital_frequency, host_mass, target_mass)
    if type(orbital_frequency) == np.ndarray:
        use_array = True

    # Get spin frequency
    if spin_frequency is None:
        if spin_period is None:
            # Assume spin-locked
            spin_frequency = orbital_frequency
        else:
            spin_frequency = days2rads(spin_period)
    if type(spin_frequency) == np.ndarray:
        use_array = True
        if type(orbital_frequency) != np.ndarray:
            orbital_frequency = orbital_frequency * np.ones_like(spin_frequency)
            semi_major_axis = semi_major_axis * np.ones_like(spin_frequency)
    elif type(orbital_frequency) == np.ndarray:
        spin_frequency = spin_frequency * np.ones_like(orbital_frequency)

    # Determine eccentricity and obliquity types
    if obliquity is None:
        obliquity = 0.
        use_obliquity = False

    if eccentricity is not None:
        if type(eccentricity) == np.ndarray:
            use_array = True
            if obliquity is not None:
                if type(obliquity) != np.ndarray:
                    obliquity = obliquity * np.ones_like(eccentricity)
            else:
                obliquity = np.zeros_like(eccentricity)
    else:
        eccentricity = 0.
        if obliquity is not None:
            if type(obliquity) == np.ndarray:
                use_array = True
                eccentricity = np.zeros_like(obliquity)

    # Get rheological model functions
    use_cpl_ctl = False
    fixed_inputs = tuple()
    if rheology.lower() in ['cpl', 'fixed_q']:
        use_cpl_ctl = True
        shear_modulus = 1.
        rheo_func = cpl_neg_imk_helper_func
        fixed_inputs = (fixed_k2, fixed_q)
    elif rheology.lower() in ['ctl']:
        use_cpl_ctl = True
        shear_modulus = 1.
        rheo_func = ctl_neg_imk_helper_func
        if fixed_dt is None:
            fixed_dt = (1. / fixed_q) * (1. / orbital_frequency)
        fixed_inputs = (fixed_k2, linear_dt, (fixed_dt,))
    else:
        # Determine strength types
        if viscosity is None or shear_modulus is None:
            raise MissingArgumentError('Viscosity and shear modulus are required for non-CPL/CTL calculations.')
        if type(viscosity) == np.ndarray:
            use_array = True
            if type(shear_modulus) != np.ndarray:
                shear_modulus = shear_modulus * np.ones_like(viscosity)
        if type(shear_modulus) == np.ndarray:
            use_array = True
            if type(viscosity) != np.ndarray:
                viscosity = viscosity * np.ones_like(viscosity)

        # Find complex compliance functions
        rheo_lookup_name = rheology.lower()
        rheo_func = known_compliance_models[rheo_lookup_name]

    # Get mode functions
    if precalculated_mode_results is None:
        tidal_term_freq_calculator, collapse_modes_func, eccentricity_func, inclination_func = \
            find_mode_manipulators(max_order_l=max_tidal_order_l,
                                   eccentricity_truncation_lvl=eccentricity_truncation_lvl,
                                   use_obliquity=use_obliquity)
    else:
        tidal_term_freq_calculator = precalculated_mode_results['tidal_term_freq_calculator']
        collapse_modes_func = precalculated_mode_results['collapse_modes_func']
        eccentricity_func = precalculated_mode_results['eccentricity_func']
        inclination_func = precalculated_mode_results['inclination_func']

    # Calculate tidal susceptibility
    tidal_susceptibility = calc_tidal_susceptibility(host_mass, target_radius, semi_major_axis)

    # Calculate eccentricity and obliquity results
    if precalculated_mode_results is not None:
        # To increase performance, you can calculate eccentricity functions outside this function
        eccentricity_results_byorderl = precalculated_mode_results['eccentricity_results']
    else:
        eccentricity_results_byorderl = eccentricity_func(eccentricity)
    obliquity_results_byorderl = inclination_func(obliquity)

    # Calculate tidal frequencies
    unique_frequencies, tidal_results_by_frequency = \
        tidal_term_freq_calculator(spin_frequency, orbital_frequency, semi_major_axis, target_radius,
                                   eccentricity_results_byorderl, obliquity_results_byorderl,
                                   multiply_modes_by_sign=True)

    # Calculate complex compliances
    if use_cpl_ctl:
        complex_compliance_by_freq = rheo_func(unique_frequencies, *fixed_inputs)
    else:
        if complex_compliance_inputs is None:
            complex_compliance_inputs = tuple()
        complex_compliance_by_freq = compliance_dict_helper(unique_frequencies,
                                                            rheo_func,
                                                            (shear_modulus**(-1), viscosity),
                                                            complex_compliance_inputs)

    # Collapse tidal modes
    tidal_heating, dUdM, dUdw, dUdO, love_number_by_orderl, negative_imk_by_orderl, effective_q_by_orderl = \
        collapse_modes_func(target_gravity, target_radius, target_density, shear_modulus, tidal_scale, host_mass,
                            tidal_susceptibility, complex_compliance_by_freq, tidal_results_by_frequency,
                            max_order_l=max_tidal_order_l, cpl_ctl_method=use_cpl_ctl)

    # Store results
    dissipation_results = {
        'tidal_heating'         : tidal_heating,
        'tidal_torque'          : host_mass * dUdO,
        'dUdM'                  : dUdM,
        'dUdw'                  : dUdw,
        'dUdO'                  : dUdO,
        'love_number_by_orderl' : love_number_by_orderl,
        'negative_imk_by_orderl': negative_imk_by_orderl,
        'effective_q_by_orderl' : effective_q_by_orderl,
        'use_array'             : use_array,
        'semi_major_axis'       : semi_major_axis,
        'orbital_frequency'     : orbital_frequency
    }

    if calculate_orbit_spin_derivatives:
        # Calculate orbit-spin derivatives
        dspin_dt = spin_rate_derivative(dUdO, target_moi, host_mass)
        da_dt, de_dt = semia_eccen_derivatives(semi_major_axis, orbital_frequency, eccentricity,
                                               target_mass, dUdM, dUdw, host_mass)
        # Store results
        dissipation_results['spin_rate_derivative'] = dspin_dt * dspin_dt_scale
        dissipation_results['eccentricity_derivative'] = de_dt * de_dt_scale
        dissipation_results['semi_major_axis_derivative'] = da_dt * da_dt_scale

    return dissipation_results


def quick_dual_body_tidal_dissipation(
        radii: Tuple[float, float], masses: Tuple[float, float],
        gravities: Tuple[float, float], densities: Tuple[float, float], mois: Tuple[float, float],
        viscosities: Union[NoneType, Tuple[NoneFloatArray, NoneFloatArray]] = None,
        shear_moduli: Union[NoneType, Tuple[NoneFloatArray, NoneFloatArray]] = None,
        rheologies: Union[str, Tuple[str, str]] = 'Maxwell',
        complex_compliance_inputs: Union[NoneType, Tuple[Tuple[float, ...], Tuple[float, ...]]] = None,
        obliquities: Union[NoneType, Tuple[NoneFloatArray, NoneFloatArray]] = None,
        spin_frequencies: Union[NoneType, Tuple[NoneFloatArray, NoneFloatArray]] = None,
        spin_periods: Union[NoneType, Tuple[NoneFloatArray, NoneFloatArray]] = None,
        tidal_scales: Tuple[float, float] = (1., 1.), fixed_k2s: Tuple[float, float] = (0.3, 0.3),
        fixed_qs: Tuple[float, float] = (100., 100.), fixed_dts: Tuple[NoneFloat, NoneFloat] = (None, None),
        eccentricity: NoneFloatArray = None,
        orbital_frequency: NoneFloatArray = None, orbital_period: NoneFloatArray = None,
        max_tidal_order_l: NoneInt = 2, eccentricity_truncation_lvl: NoneInt = 2,
        use_obliquity: Union[bool, NoneType] = True,
        da_dt_scale: float = 1., de_dt_scale: float = 1.,
        dspin_dt_scale: float = 1.) -> Dict[str, Union[FloatArray, SingleBodyResultType]]:
    """ Calculate the dual-body tidal dissipation for a target world orbiting its tidal host.

    This function pulls together several TidalPy packages to offer an easy to use interface to calculate tidal
        dissipation. However it may come with an impact to performance and flexibility. For example, there are various
        checks performed that may want to skip. It is recommended to only use this function as a starting point and
        then develop a custom function or script to suite your specific needs.

    Parameters
    ----------
    radii : Tuple[float, float]
    masses : Tuple[float, float]
    gravities : Tuple[float, float]
    densities : Tuple[float, float]
    mois : Tuple[float, float]
    viscosities : Union[NoneType, Tuple[NoneFloatArray, NoneFloatArray]] = None
    shear_moduli : Union[NoneType, Tuple[NoneFloatArray, NoneFloatArray]] = None
    rheologies : Union[str, Tuple[str, str]] = 'Maxwell'
    complex_compliance_inputs : Union[NoneType, Tuple[Tuple[float, ...], Tuple[float, ...]]] = None
    obliquities : Union[NoneType, Tuple[NoneFloatArray, NoneFloatArray]] = None
    spin_frequencies : Union[NoneType, Tuple[NoneFloatArray, NoneFloatArray]] = None
    spin_periods : Union[NoneType, Tuple[NoneFloatArray, NoneFloatArray]] = None
    tidal_scales : Tuple[float, float] = (1., 1.)
    fixed_k2s : Tuple[float, float] = (.3, .3)
    fixed_qs : Tuple[float, float] = (100., 100.)
    fixed_dts : Tuple[NoneFloat, NoneFloat] = (None, None)
        If any are `None`, a fixed_dt will be calculated as 1 / (Q * n)
    eccentricity: NoneFloatArray = None
    orbital_frequency : NoneFloatArray = None
    orbital_period : NoneFloatArray = None
    max_tidal_order_l : NoneInt = 2
    eccentricity_truncation_lvl : NoneInt = 2
    use_obliquity : Union[bool, NoneType] = True
    da_dt_scale : float = 1.
    de_dt_scale : float = 1.
    dspin_dt_scale : float = 1.

    Returns
    -------
    dissipation_results: Dict[str, Union[FloatArray, SingleBodyResultType]]
        Tidal dissipation results stored for both the `host` and `secondary`.

    """

    # Check for issues
    if type(radii) not in [list, tuple]:
        raise IncorrectArgumentType('Dual-body quick tidal dissipation calculation requires inputs in list format.')
    elif len(radii) != 2:
        raise ArgumentException('Unexpected number of list items.')

    # Flag if array functions should be used or not.
    use_array = False

    # Clean up inputs
    # Since we use the eccentricity now, outside the world calculators, we need to make sure that the obliquity and
    #    eccentricity arrays are the correct shape. The user provided use_obliquity overrides this functions best
    #    guess.
    use_obliquity_flags = list()
    obliquities_cleaned = list()
    for world_i in range(2):
        use_obliq_flag = False or use_obliquity
        clean_obliquity = 0.
        if eccentricity is not None:
            if type(eccentricity) == np.ndarray:
                clean_obliquity = np.zeros_like(eccentricity)
                use_array = True
                if obliquities is not None:
                    if obliquities[world_i] is not None:
                        use_obliq_flag = True
                        if type(obliquities[world_i]) != np.ndarray:
                            clean_obliquity = obliquities[world_i] * np.ones_like(eccentricity)
                        else:
                            clean_obliquity = obliquities[world_i]
        else:
            eccentricity = 0.
            if obliquities is not None:
                if obliquities[world_i] is not None:
                    clean_obliquity = obliquities[world_i]
                    use_obliq_flag = True
                    if type(obliquities[world_i]) == np.ndarray:
                        use_array = True
                        eccentricity = np.zeros_like(obliquities[world_i])
        obliquities_cleaned.append(clean_obliquity)
        use_obliquity_flags.append(use_obliq_flag)
    del obliquities
    use_obliquity_flags = tuple(use_obliquity_flags)
    obliquities_cleaned = tuple(obliquities_cleaned)

    # Clean up any arguments that the user did not provide
    if viscosities is None:
        viscosities = (None, None)
    if shear_moduli is None:
        shear_moduli = (None, None)
    if rheologies is None:
        rheologies = ('Maxwell', 'Maxwell')
    elif type(rheologies) is str:
        # User may have only provided one string indicating they want the same rheology for both planets.
        rheologies = (rheologies, rheologies)

    if complex_compliance_inputs is None:
        complex_compliance_inputs = (tuple(), tuple())
    elif type(complex_compliance_inputs) is tuple:
        if type(complex_compliance_inputs[0]) is tuple or type(complex_compliance_inputs[1]) is tuple:
            if complex_compliance_inputs[0] is None:
                complex_compliance_inputs = (tuple(), complex_compliance_inputs[1])
            elif complex_compliance_inputs[1] is None:
                complex_compliance_inputs = (complex_compliance_inputs[0], tuple())
        elif complex_compliance_inputs[0] is None and complex_compliance_inputs[1] is None:
            complex_compliance_inputs = (None, None)
        else:
            # User may have only provided one rheology input if they intend for the same rheology for both planets.
            complex_compliance_inputs = (complex_compliance_inputs, complex_compliance_inputs)

    if spin_frequencies is None:
        spin_frequencies = (None, None)
    if spin_periods is None:
        spin_periods = (None, None)
    if tidal_scales is None:
        tidal_scales = (1., 1.)
    if fixed_k2s is None:
        fixed_k2s = (.3, .3)
    if fixed_qs is None:
        fixed_qs = (100., 100.)
    if fixed_dts is None:
        fixed_dts = (None, None)

    # Get orbital frequency
    if orbital_frequency is None:
        if orbital_period is None:
            raise MissingArgumentError('Orbital period or frequency is required.')
        else:
            orbital_frequency = days2rads(orbital_period)
    semi_major_axis = orbital_motion2semi_a(orbital_frequency, masses[0], masses[1])
    if type(orbital_frequency) == np.ndarray:
        use_array = True

    # Determine what calculators to use (these will be the same for both the Host & Target worlds).
    precalculated_mode_results = dict()
    tidal_term_freq_calculator, collapse_modes_func, eccentricity_func, inclination_func = \
        find_mode_manipulators(max_order_l=max_tidal_order_l,
                               eccentricity_truncation_lvl=eccentricity_truncation_lvl,
                               use_obliquity=any(use_obliquity_flags))

    precalculated_mode_results['tidal_term_freq_calculator'] = tidal_term_freq_calculator
    precalculated_mode_results['collapse_modes_func'] = collapse_modes_func
    precalculated_mode_results['eccentricity_func'] = eccentricity_func
    precalculated_mode_results['inclination_func'] = inclination_func

    # Calculate the eccentricity function results as they will be the same for both the host and target worlds
    precalculated_mode_results['eccentricity_results'] = eccentricity_func(eccentricity)

    # Store dissipation results in a dict
    dissipation_results = dict()
    for world_i in range(2):
        # Assume the first world is the tidal host for naming purposes. Both world's mass will be used as the "host"
        #    mass for dissipation calculations.
        if world_i == 0:
            world_name = 'host'
            host_mass = masses[1]
        else:
            world_name = 'secondary'
            host_mass = masses[0]

        # Calculate dissipation for the world
        dissipation_results_for_world = \
            quick_tidal_dissipation(host_mass, target_radius=radii[world_i], target_mass=masses[world_i],
                                    target_gravity=gravities[world_i], target_density=densities[world_i],
                                    target_moi=mois[world_i],
                                    viscosity=viscosities[world_i], shear_modulus=shear_moduli[world_i],
                                    rheology=rheologies[world_i],
                                    complex_compliance_inputs=complex_compliance_inputs[world_i],
                                    eccentricity=eccentricity, obliquity=obliquities_cleaned[world_i],
                                    orbital_frequency=orbital_frequency, orbital_period=orbital_period,
                                    spin_frequency=spin_frequencies[world_i], spin_period=spin_periods[world_i],
                                    max_tidal_order_l=max_tidal_order_l,
                                    eccentricity_truncation_lvl=eccentricity_truncation_lvl,
                                    use_obliquity=use_obliquity_flags[world_i],
                                    tidal_scale=tidal_scales[world_i],
                                    fixed_k2=fixed_k2s[world_i], fixed_q=fixed_qs[world_i], fixed_dt=fixed_dts[world_i],
                                    precalculated_mode_results=precalculated_mode_results,
                                    calculate_orbit_spin_derivatives=False)

        # Calculate spin-rate derivatives
        dissipation_results_for_world['spin_rate_derivative'] = \
            spin_rate_derivative(dissipation_results_for_world['dUdO'], mois[world_i], host_mass) * dspin_dt_scale

        # Store Results
        dissipation_results[world_name] = dissipation_results_for_world

        # Confirm if arrays should be used or not
        if type(dissipation_results_for_world['dUdM']) == np.ndarray:
            use_array = True

    # Now that both world's dissipation is calculated we can find the dual-body evolution derivatives
    da_dt, de_dt = \
        semia_eccen_derivatives_dual(semi_major_axis, orbital_frequency, eccentricity,
                                     masses[0],
                                     dissipation_results['host']['dUdM'],
                                     dissipation_results['host']['dUdw'],
                                     masses[1],
                                     dissipation_results['secondary']['dUdM'],
                                     dissipation_results['secondary']['dUdw'])

    # Store results
    dissipation_results['eccentricity_derivative'] = de_dt * de_dt_scale
    dissipation_results['semi_major_axis_derivative'] = da_dt * da_dt_scale

    return dissipation_results


def single_dissipation_from_dict_or_world_instance(
        host: Union[dict, 'all_tidal_world_types'], secondary: Union[dict, 'all_tidal_world_types'],
        viscosity: NoneFloatArray = None, shear_modulus: NoneFloatArray = None, rheology: str = 'Maxwell',
        complex_compliance_inputs: Union[NoneType, Tuple[float, ...]] = None,
        eccentricity: NoneFloatArray = None, obliquity: NoneFloatArray = None,
        orbital_frequency: NoneFloatArray = None, orbital_period: NoneFloatArray = None,
        spin_frequency: NoneFloatArray = None, spin_period: NoneFloatArray = None,
        max_tidal_order_l: NoneInt = 2,
        eccentricity_truncation_lvl: NoneInt = 2,
        use_obliquity: Union[bool, NoneType] = True,
        tidal_scale: float = 1., fixed_k2: float = 0.3, fixed_q: float = 100.,
        fixed_dt: float = None,
        da_dt_scale: float = 1., de_dt_scale: float = 1.,
        dspin_dt_scale: float = 1.) -> SingleBodyResultType:
    """ By providing a dictionaries or TidalPy world objects, this function will pull out the necessary
        planetary parameters and calculate the single body tidal dissipation.
        It is assumed that the host world does not participate in the tidal dissipation.

    Parameters
    ----------
    host : Union[dict, 'all_tidal_world_types']
    secondary : Union[dict, 'all_tidal_world_types']
    viscosity : NoneFloatArray = None
    shear_modulus : NoneFloatArray = None
    rheology: str = 'Maxwell'
    complex_compliance_inputs : Union[NoneType, Tuple[float, ...]] = None
    eccentricity : NoneFloatArray = None
    obliquity : NoneFloatArray = None
    orbital_frequency : NoneFloatArray = None
    orbital_period : NoneFloatArray = None
    spin_frequency : NoneFloatArray = None
    spin_period : NoneFloatArray = None
    max_tidal_order_l : NoneInt = 2
    eccentricity_truncation_lvl : NoneInt = 2
    use_obliquity : Union[bool, NoneType] = True
    tidal_scale : float = 1.
    fixed_k2 : float = 0.3
    fixed_q : float = 100.
    fixed_dt : float = None
    da_dt_scale : float = 1.
    de_dt_scale : float = 1.
    dspin_dt_scale : float = 1.

    Returns
    -------
    dissipation_results : SingleBodyResultType
        Tidal dissipation results for single body dissipation.

    """

    # If the user provided a world instance, convert it to a dict.
    if type(host) == dict:
        host_dict = host
    else:
        host_dict = copy.deepcopy(host.config)
        # Pull out items that may have been calculated in the planet instance
        for param_name in ['mass']:
            if param_name in dir(host):
                host_dict[param_name] = getattr(host, param_name)
            else:
                raise MissingAttributeError
    host_mass = host_dict['mass']

    if type(secondary) == dict:
        secondary_dict = secondary
    else:
        secondary_dict = copy.deepcopy(secondary.config)
        # Pull out items that may have been calculated in the planet instance
        for param_name in ['radius', 'mass', 'gravity_surface', 'density_bulk', 'moi']:
            if param_name in dir(secondary):
                secondary_dict[param_name] = getattr(secondary, param_name)
            else:
                raise MissingAttributeError

    # Determine eccentricity and obliquity types
    if obliquity is None:
        obliquity = 0.
        use_obliquity = False

    if eccentricity is not None:
        if type(eccentricity) == np.ndarray:
            if obliquity is not None:
                if type(obliquity) != np.ndarray:
                    obliquity = obliquity * np.ones_like(eccentricity)
            else:
                obliquity = np.zeros_like(eccentricity)
    else:
        eccentricity = 0.
        if obliquity is not None:
            if type(obliquity) == np.ndarray:
                eccentricity = np.zeros_like(obliquity)

    # Calculate single body response
    single_body_response = quick_tidal_dissipation(
            host_mass=host_mass, target_radius=secondary_dict['radius'], target_mass=secondary_dict['mass'],
            target_gravity=secondary_dict['gravity_surface'], target_density=secondary_dict['density_bulk'],
            target_moi=secondary_dict['moi'],
            viscosity=viscosity, shear_modulus=shear_modulus, rheology=rheology,
            complex_compliance_inputs=complex_compliance_inputs,
            eccentricity=eccentricity, obliquity=obliquity,
            orbital_frequency=orbital_frequency, orbital_period=orbital_period,
            spin_frequency=spin_frequency, spin_period=spin_period,
            max_tidal_order_l=max_tidal_order_l,
            eccentricity_truncation_lvl=eccentricity_truncation_lvl, use_obliquity=use_obliquity,
            tidal_scale=tidal_scale, fixed_k2=fixed_k2, fixed_q=fixed_q, fixed_dt=fixed_dt,
            dspin_dt_scale=dspin_dt_scale, da_dt_scale=da_dt_scale, de_dt_scale=de_dt_scale,
            calculate_orbit_spin_derivatives=True)

    return single_body_response


def dual_dissipation_from_dict_or_world_instance(
        host: Union[dict, 'all_tidal_world_types'], secondary: Union[dict, 'all_tidal_world_types'],
        viscosities: Union[NoneType, Tuple[NoneFloatArray, NoneFloatArray]] = None,
        shear_moduli: Union[NoneType, Tuple[NoneFloatArray, NoneFloatArray]] = None,
        rheologies: Union[str, Tuple[str, str]] = 'Maxwell',
        complex_compliance_inputs: Union[NoneType, Tuple[Tuple[float, ...], Tuple[float, ...]]] = None,
        obliquities: Union[NoneType, Tuple[NoneFloatArray, NoneFloatArray]] = None,
        spin_frequencies: Union[NoneType, Tuple[NoneFloatArray, NoneFloatArray]] = None,
        spin_periods: Union[NoneType, Tuple[NoneFloatArray, NoneFloatArray]] = None,
        tidal_scales: Tuple[float, float] = (1., 1.),
        fixed_k2s: Tuple[float, float] = (0.3, 0.3),
        fixed_qs: Tuple[float, float] = (100., 100.),
        fixed_dts: Tuple[NoneFloat, NoneFloat] = (None, None),
        eccentricity: NoneFloatArray = None,
        orbital_frequency: NoneFloatArray = None,
        orbital_period: NoneFloatArray = None,
        max_tidal_order_l: NoneInt = 2, eccentricity_truncation_lvl: NoneInt = 2,
        use_obliquity: Union[bool, NoneType] = True,
        da_dt_scale: float = 1., de_dt_scale: float = 1.,
        dspin_dt_scale: float = 1.) -> Dict[str, Union[FloatArray, SingleBodyResultType]]:
    """ By providing a dictionaries or TidalPy world objects, this function will pull out the necessary
        planetary parameters and calculate the single body tidal dissipation.
        It is assumed that both the host and the secondary participate in tidal dissipation.

    Parameters
    ----------
    host : Union[dict, 'all_tidal_world_types']
    secondary : Union[dict, 'all_tidal_world_types']
    viscosities : Union[NoneType, Tuple[NoneFloatArray, NoneFloatArray]] = None
    shear_moduli : Union[NoneType, Tuple[NoneFloatArray, NoneFloatArray]] = None
    rheologies : Union[str, Tuple[str, str]] = 'Maxwell'
    complex_compliance_inputs : Union[NoneType, Tuple[Tuple[float, ...], Tuple[float, ...]]] = None
    obliquities : Union[NoneType, Tuple[NoneFloatArray, NoneFloatArray]] = None
    spin_frequencies : Union[NoneType, Tuple[NoneFloatArray, NoneFloatArray]] = None
    spin_periods : Union[NoneType, Tuple[NoneFloatArray, NoneFloatArray]] = None
    tidal_scales : Tuple[float, float] = (1., 1.)
    fixed_k2s : Tuple[float, float] = (0.3, 0.3)
    fixed_qs : Tuple[float, float] = (100., 100.)
    fixed_dts : Tuple[NoneFloat, NoneFloat] = (None, None)
    eccentricity : NoneFloatArray = None
    orbital_frequency : NoneFloatArray = None
    orbital_period : NoneFloatArray = None
    max_tidal_order_l : NoneInt = 2
    eccentricity_truncation_lvl : NoneInt = 2
    use_obliquity : Union[bool, NoneType] = True
    da_dt_scale : float = 1.
    de_dt_scale : float = 1.
    dspin_dt_scale : float = 1.

    Returns
    -------
    dissipation_results : Dict[str, Union[FloatArray, SingleBodyResultType]]
        A dictionary of the dissipation results for both the `host` and `secondary` tidal world.

    """

    # If the user provided a world instance, convert it to a dict.
    if type(host) == dict:
        host_dict = host
    else:
        host_dict = copy.deepcopy(host.config)
        # Pull out items that may have been calculated in the planet instance
        for param_name in ['radius', 'mass', 'gravity_surface', 'density_bulk', 'moi']:
            if param_name in dir(host):
                host_dict[param_name] = getattr(host, param_name)
            else:
                raise MissingAttributeError

    if type(secondary) == dict:
        secondary_dict = secondary
    else:
        secondary_dict = copy.deepcopy(secondary.config)
        # Pull out items that may have been calculated in the planet instance
        for param_name in ['radius', 'mass', 'gravity_surface', 'density_bulk', 'moi']:
            if param_name in dir(secondary):
                secondary_dict[param_name] = getattr(secondary, param_name)
            else:
                raise MissingAttributeError

    # Call the quick dual body tidal dissipation calculator
    dissipation_results = \
        quick_dual_body_tidal_dissipation(
                radii=(host_dict['radius'], secondary_dict['radius']),
                masses=(host_dict['mass'], secondary_dict['mass']),
                gravities=(host_dict['gravity_surface'], secondary_dict['gravity_surface']),
                densities=(host_dict['density_bulk'], secondary_dict['density_bulk']),
                mois=(host_dict['moi'], secondary_dict['moi']),
                viscosities=viscosities, shear_moduli=shear_moduli,
                rheologies=rheologies, complex_compliance_inputs=complex_compliance_inputs,
                obliquities=obliquities, spin_frequencies=spin_frequencies, spin_periods=spin_periods,
                tidal_scales=tidal_scales, fixed_k2s=fixed_k2s, fixed_qs=fixed_qs, fixed_dts=fixed_dts,
                eccentricity=eccentricity, orbital_frequency=orbital_frequency, orbital_period=orbital_period,
                max_tidal_order_l=max_tidal_order_l, eccentricity_truncation_lvl=eccentricity_truncation_lvl,
                use_obliquity=use_obliquity,
                da_dt_scale=da_dt_scale, de_dt_scale=de_dt_scale, dspin_dt_scale=dspin_dt_scale)

    return dissipation_results
