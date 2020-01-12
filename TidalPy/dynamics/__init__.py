from ..performance import njit

from .duel_dissipation import (eccentricity_derivative as eccentricity_derivative_duel,
                               semi_major_axis_derivative as semi_major_axis_derivative_duel)
from .single_dissipation import eccentricity_derivative, semi_major_axis_derivative, spin_rate_derivative

# Minimum tidal mode in [rad s-1]
MODE_ZERO_TOL = 1.e-12

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

from .modes_l2 import nsr_modes as nsr_modes_l2
from .modes_l2 import nsr_modes_4 as nsr_modes_l2_t4
from .modes_l2 import nsr_modes_6 as nsr_modes_l2_t6
from .modes_l2 import nsr_modes_8 as nsr_modes_l2_t8
from .modes_l2 import nsr_modes_10 as nsr_modes_l2_t10
from .modes_l2 import nsr_modes_12 as nsr_modes_l2_t12
from .modes_l2 import spin_sync_modes as spin_sync_l2
from .modes_l2 import spin_sync_modes_4 as spin_sync_l2_t4
from .modes_l2 import spin_sync_modes_6 as spin_sync_l2_t6
from .modes_l2 import spin_sync_modes_8 as spin_sync_l2_t8

from .modes_l3 import nsr_modes as nsr_modes_l3
from .modes_l3 import nsr_modes_4 as nsr_modes_l3_t4
from .modes_l3 import nsr_modes_6 as nsr_modes_l3_t6
from .modes_l3 import spin_sync_modes as spin_sync_l3
from .modes_l3 import spin_sync_modes_4 as spin_sync_l3_t4
from .modes_l3 import spin_sync_modes_6 as spin_sync_l3_t6

max_implemented_order_l = 2

mode_types = {
    # Is NSR or not
    True: {
        # order_l
        2: {
            # Truncation level
            2: nsr_modes_l2,
            4: nsr_modes_l2_t4,
            6: nsr_modes_l2_t6,
            8: nsr_modes_l2_t8,
            10: nsr_modes_l2_t10,
            12: nsr_modes_l2_t12
        },
        3: {
            # Truncation level
            2: nsr_modes_l3,
            4: nsr_modes_l3_t4,
            6: nsr_modes_l3_t6
        }
    },
    False: {
        2: {
            2: spin_sync_l2,
            4: spin_sync_l2_t4,
            6: spin_sync_l2_t6,
            8: spin_sync_l2_t8
        },
        3: {
            2: spin_sync_l3,
            4: spin_sync_l3_t4,
            6: spin_sync_l3_t6
        }
    }
}