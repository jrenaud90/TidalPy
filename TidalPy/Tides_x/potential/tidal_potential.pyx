# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""
tidal_potential.pyx
Cython/Python wrapper for TidalPy's tidal-potential truncation model hierarchy.

A tidal potential model evaluates the active tidal modes (signed frequency + the potential angular factor
U and its colatitude/longitude derivatives) at a point, for a given orbital/spin state. The 3D
stress/strain/heating kernel multiplies each mode's potential by that mode's radial response.

- SyncLowEPotential  (alias "sync_low_e") - synchronous rotation, low eccentricity, no obliquity (1 mode).
- NSRModesPotential  (alias "nsr_modes")  - moderate eccentricity, non-synchronous rotation, no obliquity
                                            (up to 9 modes n, 2n, 3n, 2o+-kn).
- NSRMedObliquityPotential (alias "nsr_modes_med_obliquity") - moderate eccentricity, moderate obliquity,
                                            non-synchronous rotation (up to 17 modes).
"""

import numpy as np
cimport numpy as cnp
cnp.import_array()

from libcpp.string cimport string
from libcpp.memory cimport unique_ptr
from libcpp.utility cimport move

from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities_x.classes_x.classes cimport PhysicsBase, c_TidalPyBaseClass

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())


cdef c_TidalPotentialConfig _build_potential_config(dict config) except *:
    """Build a c_TidalPotentialConfig from a config dict (only ``use_static`` is honored today)."""
    cdef c_TidalPotentialConfig cfg
    if config is not None and "use_static" in config:
        cfg.use_static = bool(config["use_static"])
    return cfg


cdef tuple _modeset_to_arrays(c_TidalPotentialModeSet ms):
    """Convert a C++ mode set into (modes[num], potentials[num, 6]) numpy arrays."""
    cdef int num = ms.num_modes
    cdef cnp.ndarray[cnp.float64_t, ndim=1] modes = np.empty(num, dtype=np.float64)
    cdef cnp.ndarray[cnp.float64_t, ndim=2] pots = np.empty((num, 6), dtype=np.float64)
    cdef int i
    for i in range(num):
        modes[i] = ms.mode_frequency[i]
        pots[i, 0] = ms.potential[i].U
        pots[i, 1] = ms.potential[i].dU_dtheta
        pots[i, 2] = ms.potential[i].dU_dphi
        pots[i, 3] = ms.potential[i].d2U_dtheta2
        pots[i, 4] = ms.potential[i].d2U_dphi2
        pots[i, 5] = ms.potential[i].d2U_dtheta_dphi
    return modes, pots


# =====================================================================================================================
# TidalPotentialBase
# =====================================================================================================================
cdef class TidalPotentialBase(PhysicsBase):
    """Abstract base for tidal-potential models. Instantiate a concrete subclass."""

    def __cinit__(self, *args, **kwargs):
        pass  # unique_ptr<c_TidalPotentialBase> auto-inits to nullptr

    def __init__(self, *args, **kwargs):
        raise TypeError(
            "TidalPotentialBase is abstract; instantiate a concrete model "
            "(SyncLowEPotential, NSRModesPotential).")

    def __dealloc__(self):
        self._potential_ptr.reset()
        self._ptr = NULL

    @property
    def num_modes(self) -> int:
        """Number of tidal modes this truncation produces."""
        return self._potential_ptr.get().num_modes()

    def calc_modes(
            self,
            double orbital_frequency,
            double spin_frequency,
            double eccentricity,
            double host_mass,
            double semi_major_axis,
            double radius,
            double colatitude,
            double longitude,
            double time,
            double obliquity=0.0):
        """Active tidal modes at a point.

        Returns ``(modes, potentials)``: ``modes`` is a length-``num_modes`` float64 array of signed mode
        frequencies [rad s-1]; ``potentials`` is a ``(num_modes, 6)`` float64 array whose row i holds mode
        i's ``(U, dU/dtheta, dU/dphi, d2U/dtheta2, d2U/dphi2, d2U/dtheta_dphi)``. The potential's r^2
        coefficient uses ``radius``; pass the surface radius for the 3D kernel.
        """
        cdef c_TidalPotentialState state
        state.orbital_frequency = orbital_frequency
        state.spin_frequency = spin_frequency
        state.eccentricity = eccentricity
        state.obliquity = obliquity
        state.host_mass = host_mass
        state.semi_major_axis = semi_major_axis
        cdef c_TidalPotentialModeSet ms = self._potential_ptr.get().calc_modes(
            state, radius, colatitude, longitude, time)
        return _modeset_to_arrays(ms)

    cpdef dict get_config_dict(self):
        """Return the model name as a config dict (subclasses add parameters)."""
        return {"model": self.model_name}


# =====================================================================================================================
# Tidal potential models
# =====================================================================================================================
cdef class SyncLowEPotential(TidalPotentialBase):
    """Synchronous rotation, low eccentricity, no obliquity. One mode at the orbital frequency."""

    def __cinit__(self, *args, **kwargs):
        self._sync_ptr = NULL

    def __init__(self):
        cdef c_TidalPotentialConfig config
        cdef unique_ptr[c_TidalPotentialBase] ptr = c_find_tidal_potential(c_TidalPotentialModel.SyncLowE, config)
        self._sync_ptr      = <c_SyncLowEPotential*>ptr.get()
        self._potential_ptr = move(ptr)
        self._ptr           = <c_TidalPyBaseClass*>self._potential_ptr.get()

    def __dealloc__(self):
        self._sync_ptr = NULL


cdef class NSRModesPotential(TidalPotentialBase):
    """Moderate eccentricity, non-synchronous rotation, no obliquity. Up to nine modes."""

    def __cinit__(self, *args, **kwargs):
        self._nsr_ptr = NULL

    def __init__(self, bint use_static=False):
        cdef c_TidalPotentialConfig config
        config.use_static = use_static
        cdef unique_ptr[c_TidalPotentialBase] ptr = c_find_tidal_potential(c_TidalPotentialModel.NSRModes, config)
        self._nsr_ptr       = <c_NSRModesPotential*>ptr.get()
        self._potential_ptr = move(ptr)
        self._ptr           = <c_TidalPyBaseClass*>self._potential_ptr.get()

    def __dealloc__(self):
        self._nsr_ptr = NULL

    @property
    def use_static(self) -> bool:
        """Whether the time-independent (static) potential terms are included."""
        return bool(self._nsr_ptr.get_use_static())

    cpdef dict get_config_dict(self):
        d = TidalPotentialBase.get_config_dict(self)
        d["use_static"] = bool(self._nsr_ptr.get_use_static())
        return d


cdef class NSRMedObliquityPotential(TidalPotentialBase):
    """Moderate eccentricity, moderate obliquity, non-synchronous rotation. Up to seventeen modes."""

    def __cinit__(self, *args, **kwargs):
        self._obl_ptr = NULL

    def __init__(self, bint use_static=False):
        cdef c_TidalPotentialConfig config
        config.use_static = use_static
        cdef unique_ptr[c_TidalPotentialBase] ptr = c_find_tidal_potential(
            c_TidalPotentialModel.NSRMedObliquity, config)
        self._obl_ptr       = <c_NSRMedObliquityPotential*>ptr.get()
        self._potential_ptr = move(ptr)
        self._ptr           = <c_TidalPyBaseClass*>self._potential_ptr.get()

    def __dealloc__(self):
        self._obl_ptr = NULL

    @property
    def use_static(self) -> bool:
        """Whether the time-independent (static) potential terms are included."""
        return bool(self._obl_ptr.get_use_static())

    cpdef dict get_config_dict(self):
        d = TidalPotentialBase.get_config_dict(self)
        d["use_static"] = bool(self._obl_ptr.get_use_static())
        return d


# =====================================================================================================================
# Factory
# =====================================================================================================================
def make_tidal_potential(str model_name, dict config=None) -> TidalPotentialBase:
    """Build a tidal potential model by name, returning the matching rich subclass.

    Parameters
    ----------
    model_name : str
        One of ``"sync_low_e"`` (aliases ``"synchronous_low_e"``, ``"simple"``),
        ``"nsr_modes"`` (aliases ``"nsr"``, ``"nsr_med_eccen"``), or
        ``"nsr_modes_med_obliquity"`` (aliases ``"nsr_med_obliquity"``, ``"med_obliquity"``).
    config : dict, optional
        Model parameters. Currently only ``use_static`` (bool, NSR only) is honored.

    Returns
    -------
    TidalPotentialBase
        The model subclass.

    Raises
    ------
    ValueError
        If the model name is unknown.
    """
    cdef c_TidalPotentialConfig cfg = _build_potential_config(config)
    cdef c_TidalPotentialModel model = c_tidal_potential_model_from_name(model_name.encode("utf-8"))
    cdef unique_ptr[c_TidalPotentialBase] ptr = c_find_tidal_potential(model, cfg)

    cdef SyncLowEPotential sync_potential
    cdef NSRModesPotential nsr_potential
    cdef NSRMedObliquityPotential obl_potential
    if model == c_TidalPotentialModel.SyncLowE:
        sync_potential = SyncLowEPotential.__new__(SyncLowEPotential)
        sync_potential._sync_ptr      = <c_SyncLowEPotential*>ptr.get()
        sync_potential._potential_ptr = move(ptr)
        sync_potential._ptr           = <c_TidalPyBaseClass*>sync_potential._potential_ptr.get()
        return sync_potential
    elif model == c_TidalPotentialModel.NSRModes:
        nsr_potential = NSRModesPotential.__new__(NSRModesPotential)
        nsr_potential._nsr_ptr       = <c_NSRModesPotential*>ptr.get()
        nsr_potential._potential_ptr = move(ptr)
        nsr_potential._ptr           = <c_TidalPyBaseClass*>nsr_potential._potential_ptr.get()
        return nsr_potential
    else:
        obl_potential = NSRMedObliquityPotential.__new__(NSRMedObliquityPotential)
        obl_potential._obl_ptr       = <c_NSRMedObliquityPotential*>ptr.get()
        obl_potential._potential_ptr = move(ptr)
        obl_potential._ptr           = <c_TidalPyBaseClass*>obl_potential._potential_ptr.get()
        return obl_potential
