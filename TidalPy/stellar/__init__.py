from typing import Union

from .insolation import (calc_equilibrium_temperature, equilibrium_insolation_mendez,
                         equilibrium_insolation_no_eccentricity, equilibrium_insolation_williams)

EquilibFuncType = Union[type(equilibrium_insolation_mendez), type(equilibrium_insolation_no_eccentricity),
                        type(equilibrium_insolation_williams)]

equilibrium_insolation_functions = {
    'no_eccentricity': equilibrium_insolation_no_eccentricity,
    'williams'       : equilibrium_insolation_williams,
    'mendez'         : equilibrium_insolation_mendez
    }
