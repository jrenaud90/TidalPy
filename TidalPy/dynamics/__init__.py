from .dual_dissipation import (eccentricity_derivative as eccentricity_derivative_dual,
                               semi_major_axis_derivative as semi_major_axis_derivative_dual)
from .single_dissipation import eccentricity_derivative, semi_major_axis_derivative, spin_rate_derivative

diff_eqs_dual_dissipation = {
    'eccentricity'   : eccentricity_derivative_dual,
    'semi_major_axis': semi_major_axis_derivative_dual,
    # Spin-rate is the same for both dual and single dissipation
    'spin_rate'      : spin_rate_derivative
}

diff_eqs_single_dissipation = {
    'eccentricity'   : eccentricity_derivative,
    'semi_major_axis': semi_major_axis_derivative,
    'spin_rate'      : spin_rate_derivative
}
