"""C++ global (1D) tidal dissipation models and their name-based factory.

``RheologyTide`` (alias "rheology") takes k_l from the radial solver, ``FixedQTide`` ("cpl",
"fixed_q") uses k_l*(1 - i/Q_l), ``FixedLagTide`` ("ctl", "fixed_dt") uses k_l*(1 - i*w*dt_l), and
``CTLQTide`` ("ctl_q", "fixed_dt_q") uses k_l*(1 - i*w*dt_l/Q_l). Each supplies the per-mode
dissipation multiplier -Im[k_l] used by the global mode collapse.
"""

from TidalPy.Tides_x.classes.tide import (
    TideBase,
    RheologyTide,
    FixedQTide,
    FixedLagTide,
    CTLQTide,
    make_tide,
)
from TidalPy.Tides_x.classes.collapse import collapse_global_tides

__all__ = [
    "TideBase",
    "RheologyTide",
    "FixedQTide",
    "FixedLagTide",
    "CTLQTide",
    "make_tide",
    "collapse_global_tides",
]
