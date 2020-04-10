from .duel_dissipation import (eccentricity_derivative as eccentricity_derivative_duel,
                               semi_major_axis_derivative as semi_major_axis_derivative_duel)
from .single_dissipation import eccentricity_derivative, semi_major_axis_derivative, spin_rate_derivative

diff_eqs_duel_dissipation = {
    'eccentricity'   : eccentricity_derivative_duel,
    'semi_major_axis': semi_major_axis_derivative_duel,
    # Spin-rate is the same for both duel and single dissipation
    'spin_rate'      : spin_rate_derivative
}

diff_eqs_single_dissipation = {
    'eccentricity'   : eccentricity_derivative,
    'semi_major_axis': semi_major_axis_derivative,
    'spin_rate'      : spin_rate_derivative
}
