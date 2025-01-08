""" Simple Global Approximation Tides Module
"""
from typing import Callable, Dict, TYPE_CHECKING, Tuple

from TidalPy.exceptions import (ConfigPropertyChangeError, IncorrectMethodToSetStateProperty, NotYetImplementedError,
                                UnknownModelError)
from TidalPy.utilities.performance import njit

from .base import TidesBase
from ..ctl_funcs import ctl_method_input_getters, known_ctl_methods

from TidalPy.logger import get_logger
log = get_logger("TidalPy")


if TYPE_CHECKING:
    from TidalPy.utilities.types import ComplexArray, FloatArray
    from TidalPy.structures.world_types import TidalWorldType

    from ..modes.mode_manipulation import DissipTermsArray, FreqSig, UniqueFreqType, ResultsByFreqType


# Can not cache this func since it relies on a user provided njit'd callable (ctl_method)
@njit(cacheable=False)
def ctl_neg_imk_helper_func(
    tidal_frequencies: Dict['FreqSig', 'FloatArray'], fixed_k2: float,
    ctl_method: Callable, ctl_inputs: Tuple[float, ...]
    ) -> Dict['FreqSig', 'ComplexArray']:
    """ Njit-safe helper function for calculating -Imk2 for CTL method.

    Parameters
    ----------
    tidal_frequencies : Dict[FreqSig, FloatArray]
        Njit-safe TypedDict of tidal frequencies.
    fixed_k2 : float
        World's static k2
    ctl_method : Callable
        Njit-safe CTL function
    ctl_inputs : Tuple[float, ...]
        CTL inputs

    Returns
    -------
    neg_imk_by_unique_freq : Dict[FreqSig, ComplexArray]
    """

    # Build fake values so that njit can compile the function
    fake_index = list(tidal_frequencies.keys())[0]
    fake_freq = tidal_frequencies[fake_index]
    neg_imk_by_unique_freq = {(-100, -100): fake_freq * (1. + 1.j)}

    # Real calculation
    for freq_sig, freq in tidal_frequencies.items():

        effective_q = ctl_method(freq, *ctl_inputs)
        # The 0 * fake_freq is to ensure the correct array size is used.
        neg_imk_by_unique_freq[freq_sig] = fixed_k2 * (1. - (1.j * effective_q)) + (0. * fake_freq)

        # TODO At some point numba started to not like this try and except block. Getting rid of it for now, see if an issue arises where it is needed again...
        # try:
        #     effective_q = ctl_method(freq, *ctl_inputs)
        #     # The 0 * fake_freq is to ensure the correct array size is used.
        #     neg_imk_by_unique_freq[freq_sig] = fixed_k2 * (1. - (1.j * effective_q)) + (0. * fake_freq)
        # except:
        #     # Assume that the exception was a divide by zero error due to frequency = 0.
        #     neg_imk_by_unique_freq[freq_sig] = fixed_k2 * (1. - 0.j) + (0. * fake_freq)

    # Delete fake frequency used to compile function
    del neg_imk_by_unique_freq[(-100, -100)]

    return neg_imk_by_unique_freq


@njit(cacheable=True)
def cpl_neg_imk_helper_func(tidal_frequencies: Dict['FreqSig', 'FloatArray'], fixed_k2: float, fixed_q: float) \
        -> Dict['FreqSig', 'ComplexArray']:
    """ Njit-safe helper function for calculating -Imk2 for CPL method.

    Parameters
    ----------
    tidal_frequencies : Dict[FreqSig, FloatArray]
        Njit-safe TypedDict of tidal frequencies.
    fixed_k2 : float
        World's static k2
    fixed_q : float
        Fixed dissipation Q factor.

    Returns
    -------
    neg_imk_by_unique_freq : Dict[FreqSig, ComplexArray]
    """

    # Build fake values so that njit can compile the function
    fake_index = list(tidal_frequencies.keys())[0]
    fake_freq = tidal_frequencies[fake_index]
    neg_imk_by_unique_freq = {(-100, -100): fake_freq * (1. + 1.j)}

    # Real calculation
    for freq_sig, freq in tidal_frequencies.items():
        # The 0 * fake_freq is to ensure the correct array size is used.
        neg_imk_by_unique_freq[freq_sig] = fixed_k2 * (1. - (1.j / fixed_q)) + (0. * fake_freq)

    # Delete fake frequency used to compile function
    del neg_imk_by_unique_freq[(-100, -100)]

    return neg_imk_by_unique_freq


class GlobalApproxTides(TidesBase):
    """ GlobalApproxTides
    Class used for non-layered planets (Gas Giants, Stars, Very simple homogeneous planets)

    Tides class stores model parameters and methods for heating and torque which are functions of
        (T, P, melt_frac, w, e, theata)


    Attributes
    ----------
    use_ctl
    tidal_inputs
    complex_love_by_unique_freq


    Methods
    -------
    ctl_calc_method
    ctl_calc_input_getter


    See Also
    --------
    TidalPy.tides.methods.TidesBase
    """

    model = 'global_approx'

    def __init__(self, world: 'TidalWorldType', store_config_in_world: bool = True, initialize: bool = True):
        """ Constructor for TidesBase class

        Parameters
        ----------
        world : TidalWorldType
            The world where tides are being calculated.
        store_config_in_world : bool = True
            Flag that determines if the final model's configuration dictionary should be copied into the
            `world.config` dictionary.
        initialize : bool = True
            If `True`, then an initial call to the tide's reinit method will be made at the end of construction.
        """

        super().__init__(world, store_config_in_world, initialize=False)

        # State properties
        self._tidal_inputs = None
        self._ctl_complex_love_by_unique_freq = None
        self._cpl_complex_love_by_unique_freq = None

        # Configuration properties
        self._use_ctl = None
        self._ctl_calc_method = None
        self._ctl_calc_input_getter = None

        # Call reinit for initialization
        if initialize:
            self.reinit(initial_init=True)

    def reinit(self, initial_init: bool = False):
        """ Load configurations into the Tides class and import any config-dependent functions.

        This reinit process is separate from the __init__ method because the Orbit class may need to overload some
            configurations after class initialization.

        Parameters
        ----------
        initial_init : bool = False
            This should be set to True the first time reinit is called.


        Raises
        ------
        UnknownModelError
        NotYetImplementedError

        """

        # Call parent class's reinit
        super().reinit(initial_init=initial_init)

        # Reset configuration properties
        if not initial_init:
            self._use_ctl = None
            self._ctl_calc_method = None
            self._ctl_calc_input_getter = None

        # Load in simple tide specific configurations
        self._use_ctl  = self.config['use_ctl']
        self._fixed_q  = self.config['fixed_q']
        self._fixed_k2 = self.config['static_k2']
        self._fixed_dt = self.config['fixed_dt']

        if self.use_ctl:
            # There are different CTL calculation methods that can be used.
            ctl_calc_method = self.config['ctl_calc_method']
            ctl_calc_func = None
            for method_name in [ctl_calc_method, ctl_calc_method.lower(), ctl_calc_method.title()]:
                if method_name in known_ctl_methods:
                    ctl_calc_func = known_ctl_methods[method_name]
                    ctl_calc_method = method_name
                    break

            if ctl_calc_func is None:
                raise UnknownModelError(f'Unknown CTL function requested for {self}: {ctl_calc_method}.')
            else:
                self._ctl_calc_method = ctl_calc_func

            # Build the inputs for this methods
            ctl_calc_input_signatures = ctl_method_input_getters[ctl_calc_method]

            def getter():
                _inputs = list()
                for (class_name, property_name) in ctl_calc_input_signatures:
                    if class_name in ['tides', 'self']:
                        _input = getattr(self, property_name)
                    else:
                        _input = getattr(getattr(self, class_name), property_name)
                    _inputs.append(_input)
                return tuple(_inputs)

            self._ctl_calc_input_getter = getter

        # TODO: For the simple tidal world, how to allow for higher order l? User provides k_3, k_4, ...
        if self.max_tidal_order_lvl > 2:
            raise NotYetImplementedError('Tidal order number > 2 is not implemented for global approx tides.')

    def orbit_spin_changed(
        self, eccentricity_change: bool = True, obliquity_change: bool = True,
        orbital_freq_changed: bool = True, spin_freq_changed: bool = True,
        force_obliquity_update: bool = False,
        call_world_frequency_changed: bool = True,
        call_collapse_modes: bool = True
        ) -> Tuple['UniqueFreqType', 'ResultsByFreqType']:
        """ Calculate tidal heating and potential derivative terms based on the current orbital state.

        This will also calculate new unique tidal frequencies which must then be digested by the rheological model
            at each planetary layer.


        Parameters
        ----------
        eccentricity_change : bool = True
            If there was no change in eccentricity (or if the orbit set the eccentricity) set this to False for a
            performance boost. If False, eccentricity functions won't be called.
        obliquity_change : bool = True
            If there was no change in obliquity set this to False for a performance boost.
            If False, obliquity functions won't be called.
        orbital_freq_changed : bool = True
            If there was no change in the orbital frequency set this to False for a performance boost.
        spin_freq_changed : bool = True
            If there was no change in the spin frequency set this to False.
        force_obliquity_update : bool = False
            If True, then this method will call the obliquity update function even if it would otherwise skip it
            due to obliquity dependence being turned off.
        call_world_frequency_changed : bool = True
            If True, then the method will call the world's complex compliance calculator.
            This flag is set to False for the CPL method.
        call_collapse_modes : bool = True
            If True, then this method will call collapse modes (if needed).

        Returns
        -------
        unique_tidal_frequencies : UniqueFreqType
            Each unique frequency stored as a signature (orbital motion and spin-rate coeffs),
            and the calculated frequency (combo of orbital motion and spin-rate) [rad s-1]
        tidal_terms_by_frequency : ResultsByFreqType
            Results for tidal heating, dU/dM, dU/dw, dU/dO are stored in a tuple for each tidal harmonic l and
            unique frequency.

        See Also
        --------
        TidalPy.tides.dissipation.mode_collapse

        """

        # The global approximation method does not use the planet's viscoelastic properties and does not need to ask
        #  the planet to update those properties with a new frequency. However, the user may be surprised to see
        #  the planet's complex compliances, etc. calculated for a frequency the planet is not actually in. So,
        #  we will still tell the planet to make its updates even though they are not needed for this Tides class.
        call_world_frequency_changed = call_world_frequency_changed

        # Find all the tidal modes using the base class' method. Prevent it from calling the world's complex compliance
        #   calculator.
        super().orbit_spin_changed(
            eccentricity_change=eccentricity_change, obliquity_change=obliquity_change,
            orbital_freq_changed=orbital_freq_changed, spin_freq_changed=spin_freq_changed,
            force_obliquity_update=force_obliquity_update,
            call_world_frequency_changed=call_world_frequency_changed, call_collapse_modes=False
            )

        # If the CTL method is used then the dissipation efficiency will change with frequency.
        #    In CPL: Dissipation ~ k_2 / Q
        #    In CTL: Dissipation ~ k_2 / Q * (1 / \Delta{}t w) (see Correia 2009 and Heller+2011)
        #    \Delta{}t is often set equal to 1 so that CTL dissipation ~ k_2 / (Q * w)
        #    w is an ill-defined frequency. Generally it is set to the orbital motion, but some set it to the spin-rate
        #        for a world experiencing NSR (see Correia 2009).
        if self._new_tidal_frequencies:
            if self.use_ctl:
                # CTL Method
                # Get CTL inputs
                # OPT: These getters could be replaced by a set_fixed_q or set_fixed_dt since they really won't change
                #   often. It is a waste of resources to keep calling these getters.
                ctl_inputs = self.ctl_calc_input_getter()

                # Calculate new values
                self._ctl_complex_love_by_unique_freq = \
                    ctl_neg_imk_helper_func(
                        self.unique_tidal_frequencies, self.fixed_k2,
                        self.ctl_calc_method, ctl_inputs
                        )
            else:
                # CPL Method
                self._cpl_complex_love_by_unique_freq = \
                    cpl_neg_imk_helper_func(
                        self.unique_tidal_frequencies, self.fixed_k2,
                        self.fixed_q
                        )

        if self._need_to_collapse_modes and call_collapse_modes:
            self.collapse_modes()

        # Return frequencies and tidal terms
        return self.unique_tidal_frequencies, self.tidal_terms_by_frequency

    def fixed_q_dt_changed(self):
        """ The fixed tidal dissipation parameters (fixed-q or fixed-dt) have changed. Make any necessary updates. """

        super().fixed_q_dt_changed()

        self.collapse_modes()

    def clear_state(self):
        """ Clear the state for the global tides model """

        super().clear_state()

        self._tidal_inputs = None
        self._ctl_complex_love_by_unique_freq = None
        self._cpl_complex_love_by_unique_freq = None

    def collapse_modes(self) -> 'DissipTermsArray':
        """ Calculate Global Love number based on current thermal state.

        Requires a prior orbit_spin_changed() call as unique frequencies are used to calculate the complex compliances
            used to calculate the Love numbers.

        Returns
        -------
        tidal_heating : np.ndarray
            Tidal heating [W]
            This could potentially restricted to a layer or for an entire planet.
        dUdM : np.ndarray
            Tidal potential derivative with respect to the mean anomaly [J kg-1 radians-1]
            This could potentially restricted to a layer or for an entire planet.
        dUdw : np.ndarray
            Tidal potential derivative with respect to the pericentre [J kg-1 radians-1]
            This could potentially restricted to a layer or for an entire planet.
        dUdO : np.ndarray
            Tidal potential derivative with respect to the planet's node [J kg-1 radians-1]
            This could potentially restricted to a layer or for an entire planet.

        See Also
        --------
        TidalPy.tides.Tides.orbit_spin_changed
        """

        super().collapse_modes()

        # Pull out parameters that are used multiple times in this method
        tidal_terms_by_frequency = self.tidal_terms_by_frequency
        complex_love_by_unique_freq = self.complex_love_by_unique_freq

        # Check to see if all the needed state properties are present and then begin calculations
        if tidal_terms_by_frequency is not None and complex_love_by_unique_freq is not None:
            tidal_scale, radius, bulk_density, gravity_surf = self.tidal_inputs

            # Shear modulus is not used in the CTL/CPL scheme. However, it needs to be provided as a number to the
            #    collapse_modes function.
            shear_modulus = 1.

            # Mode collapse will parse through tidal order-l and all unique frequencies and calculate global dissipation
            #    values
            tidal_heating, dUdM, dUdw, dUdO, love_number_by_orderl, negative_imk_by_orderl, effective_q_by_orderl = \
                self.collapse_modes_func(
                    gravity_surf, radius, bulk_density, shear_modulus,
                    tidal_scale,
                    self.tidal_host.mass, self.tidal_susceptibility,
                    complex_love_by_unique_freq,
                    tidal_terms_by_frequency, self.max_tidal_order_lvl,
                    cpl_ctl_method=True
                    )

            # Calculation finished. Store info in accessible containers
            self._tidal_heating_global = tidal_heating
            self._dUdM = dUdM
            self._dUdw = dUdw
            self._dUdO = dUdO
            self._global_love_by_orderl = love_number_by_orderl
            self._global_negative_imk_by_orderl = negative_imk_by_orderl
            self._effective_q_by_orderl = effective_q_by_orderl

            # Tell the parent class to update the world's dissipation flag.
            self.world.dissipation_changed()

        # Return tidal heating and derivatives
        return self.tidal_heating_global, self.dUdM, self.dUdw, self.dUdO

    # # Configuration properties
    @property
    def use_ctl(self) -> bool:
        """ Flag set to `True` if the CTL method is being used (over the CPL method) """
        return self._use_ctl

    @use_ctl.setter
    def use_ctl(self, value):
        raise ConfigPropertyChangeError

    @property
    def ctl_calc_method(self) -> Callable:
        """ The method used to calculate the CTL method's frequency dependence """
        return self._ctl_calc_method

    @ctl_calc_method.setter
    def ctl_calc_method(self, value):
        raise ConfigPropertyChangeError

    @property
    def ctl_calc_input_getter(self) -> Callable:
        """ Functions used to find the inputs for the ctl_calc_method """
        return self._ctl_calc_input_getter

    @ctl_calc_input_getter.setter
    def ctl_calc_input_getter(self, value):
        raise ConfigPropertyChangeError

    # # State properties
    @property
    def tidal_inputs(self) -> Tuple[float, float, float, float]:
        """ The inputs required to calculate tides - these could change dynamically so they need to be pulled live """
        return self.world.tidal_scale, self.radius, self.world.density_bulk, self.world.gravity_surface

    @tidal_inputs.setter
    def tidal_inputs(self, value):
        raise IncorrectMethodToSetStateProperty

    @property
    def complex_love_by_unique_freq(self):
        """ -Im[k2] stored by unique frequency signature. Chooses between CPL and CTL """
        if self.use_ctl:
            return self._ctl_complex_love_by_unique_freq
        else:
            return self._cpl_complex_love_by_unique_freq

    @complex_love_by_unique_freq.setter
    def complex_love_by_unique_freq(self, value):
        raise IncorrectMethodToSetStateProperty
