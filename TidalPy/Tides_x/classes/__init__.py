"""TidalPy Tides_x.classes — C++ global (1D) tidal dissipation model hierarchy.

Exposes the four tide models and a name-based factory:

- ``RheologyTide``  (alias ``"rheology"``)               — k_l from the radial solver.
- ``FixedQTide``    (alias ``"cpl"``/``"fixed_q"``)      — constant phase lag, k_l*(1 - i/Q_l).
- ``FixedLagTide``  (alias ``"ctl"``/``"fixed_dt"``)     — constant time lag, k_l*(1 - i*w*dt_l).
- ``CTLQTide``      (alias ``"ctl_q"``/``"fixed_dt_q"``) — k_l*(1 - i*w*dt_l/Q_l).

Each model maps a per-mode Love number to the dissipation multiplier ``-Im[k_l]`` used by
the global mode collapse (``TidalPy.Tides_x.potential.collapse_global_tides``).
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
