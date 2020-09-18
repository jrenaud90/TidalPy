""" Simple Tides Module
"""

from typing import TYPE_CHECKING, Dict, Tuple

import numpy as np

from .base import TidesBase
from .defaults import tide_defaults
from ..eccentricityFuncs import EccenOutput
from ..mode_manipulation import FreqSig, DissipTermsArray
from ...exceptions import (NotYetImplementedError, ParameterValueError,
                           OuterscopePropertySetError, FailedForcedStateUpdate, IncorrectMethodToSetStateProperty)
from ...utilities.types import FloatArray

if TYPE_CHECKING:
    from ...structures.worlds import TidalWorldType


class SimpleTides(TidesBase):
    """ SimpleTides
    Class used for non-layered planets (Gas Giants, Stars, Very simple homogeneous planets)

    Tides class stores model parameters and methods for heating and torque which are functions of
        (T, P, melt_frac, w, e, theata)

    See Also
    --------
    TidalPy.tides.classes.TidesBase
    """

    model = 'simple'
    default_config = tide_defaults['simple']

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

        super().__init__(world, store_config_in_world, initialize=initialize)

        # Ensure the tidal order and orbital truncation levels make sense
        # TODO: For the simple tidal world, how to allow for higher order l? User provides k_3, k_4, ...
        if self.max_tidal_order_lvl > 2:
            raise NotYetImplementedError(f'Tidal order {self.max_tidal_order_lvl} has not been implemented for '
                                            'simple tidal worlds yet.')
        if self.eccentricity_truncation_lvl % 2 != 0:
            raise ParameterValueError('Orbital truncation level must be an even integer.')
        if self.eccentricity_truncation_lvl <= 2:
            raise ParameterValueError('Orbital truncation level must be greater than or equal to 2.')
        if self.eccentricity_truncation_lvl not in (2, 4, 6, 8, 10, 12, 14, 16, 18, 20):
            if self.max_tidal_order_lvl == 2 and self.eccentricity_truncation_lvl == 22:
                # This was implemented in v0.2.1
                pass
            else:
                raise NotYetImplementedError(f'Orbital truncation level of {self.eccentricity_truncation_lvl} has not '
                                              f'been implemented yet.')

        # State properties
        self._tidal_inputs = None
        self._neg_imk_ctl_by_unique_freq = None
        self._neg_imk_cpl_by_unique_freq = None

        # Flags
        self.use_ctl = self.config['use_ctl']
        self._thermal_set = True  # Simple Tide class does not care about thermal state

    def clear_state(self):
        """ Clear the state for the simple tides model """

        super().clear_state()

        self._neg_imk_ctl_by_unique_freq = None
        self._neg_imk_cpl_by_unique_freq = None

    def update_orbit_spin(self, force_update: bool = True, eccentricity_change: bool = True,
                          obliquity_change: bool = True, frequency_change: bool = True,
                          eccentricity_results: Dict[int, EccenOutput] = None) -> \
            Tuple[Dict[FreqSig, FloatArray], Dict[FreqSig, Dict[int, DissipTermsArray]]]:
        """ Calculate tidal heating and potential derivative terms based on the current orbital state.

        This will also calculate new unique tidal frequencies which must then be digested by the rheological model
            at each planetary layer.


        Parameters
        ----------
        force_update : bool = True
            If True, will raise an error if the update does not successful complete.
            Failed completions are usually a result of a missing state property(ies).
        eccentricity_change : bool = True
            If there was no change in eccentricity (or if the orbit set the eccentricity) set this to False for a
            performance boost. If False, eccentricity functions won't be called.
        obliquity_change : bool = True
            If there was no change in obliquity set this to False for a performance boost.
            If False, obliquity functions won't be called.
        frequency_change : bool = True
            If there was no change in orbital or rotation frequency set this to False for a performance boost.
            If False, calculate_terms won't be called.
        eccentricity_results : Dict[int, EccenOutput] = None
            Eccentricity functions can be calculated by the Orbit class (or by the user). Pass a pre-caclualted result
            for a performance boost.

        Returns
        -------
        unique_tidal_frequencies : Dict[FreqSig, np.ndarray]
            Each unique frequency stored as a signature (orbital motion and spin-rate coeffs), and the calculated frequency
                (combo of orbital motion and spin-rate) [rad s-1]
        tidal_terms_by_frequency : Dict[FreqSig, Dict[int, DissipTermsArray]]
            Results for tidal heating, dU/dM, dU/dw, dU/dO are stored in a tuple for each tidal harmonic l and
                unique frequency.

        See Also
        --------
        TidalPy.tides.dissipation.mode_collapse
        """

        # If the CTL method is used then the dissipation efficiency will change with frequency.
        #    In CPL: Dissipation ~ k_2 / Q
        #    In CTL: Dissipation ~ k_2 / Q * (1 / \Delta{}t w) (see Correia 2009 and Heller+2011)
        #    \Delta{}t is often set equal to 1 so that CTL dissipation ~ k_2 / (Q * w)
        #    w is a ill-defined frequency. Generally it is set to the orbital motion, but some set it to the spin-rate
        #        for a world experiencing NSR (see Correia 2009).

        if self.unique_tidal_frequencies is not None:

            real_val = self.fixed_k2

            if self.use_ctl:
                # Calculate new values
                self._neg_imk_ctl_by_unique_freq = \
                    {freq_sig: real_val + 1.0j * self.fixed_k2 / (self.fixed_q * freq * self.fixed_dt)
                     for freq_sig, freq in self.unique_tidal_frequencies.items()}
            else:

                self._neg_imk_cpl_by_unique_freq = \
                    {freq_sig: real_val + 1.0j * self.fixed_k2 / self.fixed_q
                     for freq_sig, freq in self.unique_tidal_frequencies.items()}

        # Return frequencies and tidal terms
        return self.unique_tidal_frequencies, self.tidal_terms_by_frequency

    def collapse_modes(self, force_update: bool = True) -> DissipTermsArray:
        """ Calculate Global Love number based on current thermal state.

        Requires a prior update_orbit_spin() call as unique frequencies are used to calculate the complex compliances
            used to calculate the Love numbers.

        Parameters
        ----------
        force_update : bool = True
            If True, will raise an error if the update does not successful complete.
            Failed completions are usually a result of a missing state property(ies).

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
        TidalPy.tides.Tides.update_orbit_spin
        """

        # Check to see if all the needed state properties are present
        all_values_present = True

        # Check to see if tidal terms have been loaded.
        if self.tidal_terms_by_frequency is None:
            # Attempt to update the orbit to see if we can load them
            self.update_orbit_spin(force_update=False)

            if self.tidal_terms_by_frequency is None:
                # Update did not help.
                all_values_present = False

        if all_values_present:
            tidal_scale, radius, bulk_density, gravity_surf = self.tidal_inputs

            # Shear modulus is not used in the CTL/CPL scheme. However, it needs to be provided as a number to the
            #    collapse_modes function.
            shear_modulus = 1.

            # Mode collapse will parse through tidal order-l and all unique frequencies and calculate global dissipation
            #    values
            tidal_heating, dUdM, dUdw, dUdO, love_number, negative_imk = \
                self.collapse_modes_func(gravity_surf, radius, bulk_density, shear_modulus,
                                         self.neg_imk_by_unique_freq,
                                         self.tidal_terms_by_frequency, self.tidal_susceptibility,
                                         self.tidal_host.mass,
                                         tidal_scale, cpl_ctl_method=True)

            # Calculation finished. Store info in accessible containers
            self._tidal_heating_global = tidal_heating
            self._dUdM = dUdM
            self._dUdw = dUdw
            self._dUdO = dUdO
            self._negative_imk_global = negative_imk

            # Now tell other methods to update now that derivatives and heating has been altered
            # TODO: orbit derivatives
            self.set_spin_derivative()

            # Return tidal heating and derivatives
            return self.tidal_heating_global, self.dUdM, self.dUdw, self.dUdO

        else:
            if force_update:
                raise FailedForcedStateUpdate


    # # State properties
    @property
    def tidal_inputs(self) -> Tuple[float, float, float, float]:
        """ The inputs required to calculate tides - these could change dynamically so they need to be pulled live """
        return self.world.tidal_scale, self.radius, self.world.density_bulk, self.world.gravity_surface

    @tidal_inputs.setter
    def tidal_inputs(self, value):
        raise IncorrectMethodToSetStateProperty

    @property
    def neg_imk_ctl_by_unique_freq(self) -> Dict[FreqSig, np.ndarray]:
        """ -Im[k2] stored by unique frequency signature. Used for the CTL method """
        return self._neg_imk_ctl_by_unique_freq

    @neg_imk_ctl_by_unique_freq.setter
    def neg_imk_ctl_by_unique_freq(self, value):
        raise IncorrectMethodToSetStateProperty

    @property
    def neg_imk_cpl_by_unique_freq(self) -> Dict[FreqSig, np.ndarray]:
        """ -Im[k2] stored by unique frequency signature. Used for the CPL method """
        return self._neg_imk_cpl_by_unique_freq

    @neg_imk_cpl_by_unique_freq.setter
    def neg_imk_cpl_by_unique_freq(self, value):
        raise IncorrectMethodToSetStateProperty

    @property
    def neg_imk_by_unique_freq(self):
        """ -Im[k2] stored by unique frequency signature. Chooses between CPL and CTL """
        if self.use_ctl:
            return self.neg_imk_ctl_by_unique_freq
        else:
            return self.neg_imk_cpl_by_unique_freq

    @neg_imk_by_unique_freq.setter
    def neg_imk_by_unique_freq(self, value):
        raise IncorrectMethodToSetStateProperty


    # # Outer-scope properties
    @property
    def fixed_k2(self) -> float:
        """ Outer-scope wrapper for world.static_love """
        return self.world.static_love

    @fixed_k2.setter
    def fixed_k2(self, value):
        raise OuterscopePropertySetError

    @property
    def fixed_q(self) -> float:
        """ Outer-scope wrapper for world.fixed_q """
        return self.world.fixed_q

    @fixed_q.setter
    def fixed_q(self, value):
        raise OuterscopePropertySetError

    @property
    def fixed_dt(self):
        """ Outer-scope wrapper for world.fixed_dt """
        return self.world.fixed_dt

    @fixed_dt.setter
    def fixed_dt(self, value):
        raise OuterscopePropertySetError
