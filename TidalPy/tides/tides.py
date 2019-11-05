import copy
from functools import partial
from typing import Tuple

import numpy as np
from . import calculate_tides
from ..exceptions import (AttributeNotSetError, ImplementationException, ParameterMissingError, UnknownModelError,
                          BadValueError)
from ..rheology import andrade_frequency_models, compliance_models
from ..rheology.compliance import ComplianceModelSearcher
from ..rheology.defaults import rheology_param_defaults
from ..types import FloatArray, ArrayNone
from ..utilities.model import LayerModel


FAKE_FREQS = (np.asarray([0.]),)


class Tides(LayerModel):
    default_config = copy.deepcopy(rheology_param_defaults)
    config_key = 'rheology'

    def __init__(self, layer):

        super().__init__(layer=layer, function_searcher=None, call_reinit=True)

        # Override model if layer has tidal off
        if not self.layer.config['is_tidal']:
            self.model = 'off'

        # Pull out model information
        self.orbital_truncation_level = self.world.config['orbital_truncation_level']
        self.order_l = self.world.config['tidal_order_l']

        # Pull out information about the planet/layer
        self.config['planet_beta'] = self.layer.world.beta
        self.config['quality_factor'] = self.layer.world.fixed_q

        if int(self.order_l) != self.order_l or self.order_l < 2:
            raise BadValueError('Tidal order l must be an integer >= 2')
        if self.order_l > 3:
            raise ImplementationException('Tidal order l > 3 has not been fully implemented in TidalPy.')

        if int(self.orbital_truncation_level) != self.orbital_truncation_level or \
                self.orbital_truncation_level % 2 != 0 or self.orbital_truncation_level < 2:
            raise BadValueError('Orbital truncation level must be an even integer >= 2')
        if self.orbital_truncation_level > 6:
            raise ImplementationException('Orbital truncation level > 6 has not been fully implemented in TidalPy.')

        # Setup Love number calculator
        self.full_calculation = True
        if self.model == 'regular':
            # FIXME: I don't like how the model is stored for the tidal calculations.
            #    right now this stuff is stored in the /rheology/default.py which doesn't make a lot of sense.
            pass
        elif self.model == 'off':
            self.full_calculation = False
        elif self.model == 'multilayer':
            # TODO: Multilayer code
            raise ImplementationException
        else:
            raise UnknownModelError

        # Setup complex compliance calculator
        model_searcher = ComplianceModelSearcher(compliance_models, andrade_frequency_models, self.config,
                                                 defaults_require_key=False)

        if not self.full_calculation:
            self.rheology_name = 'off'
        else:
            self.rheology_name = self.config['compliance_model']

        comp_func, comp_input, comp_live_args = model_searcher.find_model(self.rheology_name)
        # TODO: As of 0.1.0a no compliance functions have live args, so comp_live_args is never used.
        self.compliance_func = comp_func
        self.compliance_inputs = comp_input

        # Switches
        self.is_spin_sync = self.layer.is_spin_sync

    def _calculate(self) -> Tuple[ArrayNone, ArrayNone]:
        """ Calculates the Tidal Heating and Tidal Torque for this layer.

        Returns
        -------
        total_heating : ArrayNone
            Total Tidal Heating [W]
        total_torque : ArrayNone
            Total Tidal Torque [N m]

        """

        # Physical and Frequency Data
        try:
            shear = self.layer.shear_modulus
            viscosity = self.layer.viscosity
            orbital_freq = self.world.orbital_freq
            spin_freq = self.world.spin_freq
            eccentricity = self.world.eccentricity
            inclination = self.world.inclination
        except AttributeNotSetError:
            raise ParameterMissingError
        if shear is None:
            raise ParameterMissingError
        if viscosity is None:
            raise ParameterMissingError
        if orbital_freq is None:
            raise ParameterMissingError
        if spin_freq is None:
            if self.is_spin_sync:
                spin_freq = orbital_freq
            else:
                raise ParameterMissingError
        if eccentricity is None:
            eccentricity = np.asarray([0])
        if inclination is None:
            inclination = np.asarray([0])

        if self.full_calculation:
            # Make a call to the calculate tides function
            tidal_heating, tidal_torque = \
                calculate_tides(shear, viscosity, self.layer.gravity, self.layer.radius, self.layer.density,
                                self.world.tidal_susceptibility, self.compliance_func, self.compliance_inputs,
                                self.world.semi_major_axis, eccentricity, inclination,
                                orbital_freq, spin_freq,
                                tidal_volume_fraction=self.layer.tidal_scale, use_nsr=self.is_spin_sync,
                                truncation=self.orbital_truncation_level, order_l=self.order_l)
        else:
            tidal_heating, tidal_torque = None, None

        return tidal_heating, tidal_torque
