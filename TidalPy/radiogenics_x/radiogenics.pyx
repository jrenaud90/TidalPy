# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""
radiogenics.pyx
Cython/Python wrappers for TidalPy's radiogenics model hierarchy.

Exposes the three radiogenic-heating models:

- ``OffRadiogenics``     (alias ``"none"``)     — radiogenics disabled (heating == 0).
- ``IsotopeRadiogenics``                        — sum of individual decaying isotopes.
- ``FixedRadiogenics``   (alias ``"constant"``) — single lumped rate with optional decay.

Each model computes the total radiogenic heating [W] via ``calc_heating(time_s,
mass_kg)``. All quantities are MKS: time and half-lives in seconds, mass in kg,
heat production rates in W/kg, heating in W.

References
----------
- Hussmann and Spohn (2004); Turcotte and Schubert (2001) — chondritic isotope data.
- Castillo-Rogez et al. (2007) — long- and short-lived radiogenic isotopes.
"""

from libcpp cimport bool as cpp_bool
from libcpp.memory cimport unique_ptr
from libcpp.string cimport string
from libcpp.utility cimport move
from libcpp.vector cimport vector

import numpy as np

from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.Utilities_x.classes_x.classes cimport PhysicsBase, c_TidalPyBaseClass

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())


# Seconds in one mega-year (Julian year, 365.25 days). Used to convert the
# global-config isotope data (stored in Myr) to MKS seconds.
_SECONDS_PER_MYR = 1.0e6 * 365.25 * 24.0 * 3600.0


# =====================================================================================================================
# Internal helpers for vectorized solving
# =====================================================================================================================
cdef void _fill_vector(double[::1] src, vector[double]& dst) noexcept:
    """Copy a contiguous 1-D float64 memoryview into a std::vector[double]."""
    cdef Py_ssize_t n = src.shape[0]
    cdef Py_ssize_t i
    dst.resize(n)
    for i in range(n):
        dst[i] = src[i]


cdef object _double_vector_to_ndarray(vector[double]& src, tuple shape):
    """Build a float64 ndarray (of the given shape) from a std::vector."""
    cdef Py_ssize_t n = <Py_ssize_t>src.size()
    cdef Py_ssize_t i
    out = np.empty(n, dtype=np.float64)
    cdef double[::1] mv = out
    for i in range(n):
        mv[i] = src[i]
    return out.reshape(shape)


cdef void _build_isotopes(object heat_production_w_kg, object half_lives_s,
                          object mass_fracs, object concentrations,
                          object names, vector[c_Isotope]& dst):
    """Build a std::vector[c_Isotope] from parallel Python sequences.

    The four numeric arrays must have the same length; ``names`` is an optional
    sequence of isotope labels (auto-generated as ``isotope_<i>`` when ``None``).
    """
    cdef double[::1] a_hpr  = np.ascontiguousarray(heat_production_w_kg, dtype=np.float64).ravel()
    cdef double[::1] a_half = np.ascontiguousarray(half_lives_s,         dtype=np.float64).ravel()
    cdef double[::1] a_frac = np.ascontiguousarray(mass_fracs,           dtype=np.float64).ravel()
    cdef double[::1] a_conc = np.ascontiguousarray(concentrations,       dtype=np.float64).ravel()
    cdef Py_ssize_t n = a_hpr.shape[0]
    if a_half.shape[0] != n or a_frac.shape[0] != n or a_conc.shape[0] != n:
        raise ValueError(
            "TidalPy: isotope heat_production, half_lives, mass_fracs, and "
            "concentrations must all have the same length.")
    if names is not None and len(names) != n:
        raise ValueError("TidalPy: isotope 'names' length must match the isotope arrays.")
    dst.clear()
    cdef Py_ssize_t i
    cdef string nm
    for i in range(n):
        if names is not None:
            nm = str(names[i]).encode("utf-8")
        else:
            nm = ("isotope_%d" % i).encode("utf-8")
        dst.push_back(c_Isotope(nm, a_hpr[i], a_half[i], a_frac[i], a_conc[i]))


cdef object _isotopes_to_arrays(const vector[c_Isotope]& isotopes):
    """Extract (heat_production, half_lives, mass_fracs, concentrations, names).

    The four numeric outputs are float64 ndarrays; ``names`` is a list of str.
    """
    cdef Py_ssize_t n = <Py_ssize_t>isotopes.size()
    cdef Py_ssize_t i
    hpr  = np.empty(n, dtype=np.float64)
    half = np.empty(n, dtype=np.float64)
    frac = np.empty(n, dtype=np.float64)
    conc = np.empty(n, dtype=np.float64)
    cdef double[::1] m_hpr  = hpr
    cdef double[::1] m_half = half
    cdef double[::1] m_frac = frac
    cdef double[::1] m_conc = conc
    names = []
    for i in range(n):
        m_hpr[i]  = isotopes[i].heat_production_w_kg
        m_half[i] = isotopes[i].half_life_s
        m_frac[i] = isotopes[i].mass_frac
        m_conc[i] = isotopes[i].concentration
        names.append(isotopes[i].name.decode("utf-8"))
    return hpr, half, frac, conc, names


cdef object _solve_heating(c_RadiogenicsBase* model, object time, object mass):
    """Solve radiogenic heating for float and/or np.ndarray inputs.

    Picks the most specific C++ vectorized routine for the input pattern:

    - both scalar              -> scalar calc_heating (returns float)
    - time array, mass scalar  -> vectorize_time
    - mass array, time scalar  -> vectorize_mass
    - any other mix of arrays  -> broadcast both -> vectorize_all

    Returns a Python ``float`` for all-scalar input, else an ``np.ndarray``
    (float64) broadcast to the common shape.
    """
    cdef cpp_bool t_arr = isinstance(time, np.ndarray)
    cdef cpp_bool m_arr = isinstance(mass, np.ndarray)

    cdef vector[double] vtime, vmass
    cdef vector[double] vout
    cdef double[::1] mv

    # All scalar -> scalar result.
    if not (t_arr or m_arr):
        return float(model.calc_heating(<double>time, <double>mass))

    # Time varies; mass constant.
    if t_arr and not m_arr:
        time_arr = np.ascontiguousarray(time, dtype=np.float64)
        mv = time_arr.ravel()
        _fill_vector(mv, vtime)
        model.calc_heating_vectorize_time(vtime, <double>mass, vout)
        return _double_vector_to_ndarray(vout, time_arr.shape)

    # Mass varies; time constant.
    if m_arr and not t_arr:
        mass_arr = np.ascontiguousarray(mass, dtype=np.float64)
        mv = mass_arr.ravel()
        _fill_vector(mv, vmass)
        model.calc_heating_vectorize_mass(<double>time, vmass, vout)
        return _double_vector_to_ndarray(vout, mass_arr.shape)

    # General case: broadcast both and vary everything.
    t_b, m_b = np.broadcast_arrays(
        np.asarray(time, dtype=np.float64),
        np.asarray(mass, dtype=np.float64))
    t_c = np.ascontiguousarray(t_b)
    m_c = np.ascontiguousarray(m_b)
    mv = t_c.ravel(); _fill_vector(mv, vtime)
    mv = m_c.ravel(); _fill_vector(mv, vmass)
    model.calc_heating_vectorize_all(vtime, vmass, vout)
    return _double_vector_to_ndarray(vout, t_c.shape)


# =====================================================================================================================
# RadiogenicsBase
# =====================================================================================================================
cdef class RadiogenicsBase(PhysicsBase):
    """Abstract base for all radiogenics models.

    Holds the owning ``unique_ptr`` to the most-derived C++ model object and
    exposes the shared radiogenic-heating calculations. Not directly
    instantiable — construct a concrete model (e.g. ``IsotopeRadiogenics``).
    """

    def __cinit__(self, *args, **kwargs):
        pass  # unique_ptr<c_RadiogenicsBase> auto-inits to nullptr; concrete models set it

    def __init__(self, *args, **kwargs):
        raise TypeError(
            "RadiogenicsBase is abstract; instantiate a concrete model "
            "(OffRadiogenics, IsotopeRadiogenics, FixedRadiogenics)."
        )

    def __dealloc__(self):
        # unique_ptr frees the most-derived C++ object; null the base observer ptr.
        self._radiogenics_ptr.reset()
        self._ptr = NULL

    # ------------------------------------------------------------------------------------------------------------------
    # Calculations
    # ------------------------------------------------------------------------------------------------------------------
    def calc_heating(self, double time_s, double mass_kg) -> float:
        """Total radiogenic heating [W] for the given time and mass.

        Parameters
        ----------
        time_s : float
            Elapsed time [s] (shares its zero point with the model reference time).
        mass_kg : float
            Mass of the radiogenic material [kg].

        Returns
        -------
        float
            Radiogenic heating [W].

        Assumptions
        -----------
        - Exponential radioactive decay from a reference time.
        - All inputs and outputs are MKS.
        """
        self._check_ptr()
        return self._radiogenics_ptr.get().calc_heating(time_s, mass_kg)

    def calc_heating_vectorize_time(self, time_s, double mass_kg):
        """Radiogenic heating over a time sweep at constant mass.

        Parameters
        ----------
        time_s : array_like
            1-D sequence of elapsed times [s].
        mass_kg : float
            Constant radiogenic mass [kg].

        Returns
        -------
        numpy.ndarray
            Radiogenic heating [W] for each time (float64, same length).
        """
        self._check_ptr()
        cdef vector[double] vtime
        cdef vector[double] vout
        cdef double[::1] mv
        time_c = np.ascontiguousarray(time_s, dtype=np.float64).ravel()
        mv = time_c; _fill_vector(mv, vtime)
        self._radiogenics_ptr.get().calc_heating_vectorize_time(vtime, mass_kg, vout)
        return _double_vector_to_ndarray(vout, time_c.shape)

    def calc_heating_vectorize_mass(self, double time_s, mass_kg):
        """Radiogenic heating over a mass sweep at constant time.

        Parameters
        ----------
        time_s : float
            Constant elapsed time [s].
        mass_kg : array_like
            1-D sequence of radiogenic masses [kg].

        Returns
        -------
        numpy.ndarray
            Radiogenic heating [W] for each mass (float64, same length).
        """
        self._check_ptr()
        cdef vector[double] vmass
        cdef vector[double] vout
        cdef double[::1] mv
        mass_c = np.ascontiguousarray(mass_kg, dtype=np.float64).ravel()
        mv = mass_c; _fill_vector(mv, vmass)
        self._radiogenics_ptr.get().calc_heating_vectorize_mass(time_s, vmass, vout)
        return _double_vector_to_ndarray(vout, mass_c.shape)

    def calc_heating_vectorize_all(self, time_s, mass_kg):
        """Radiogenic heating over element-wise (time, mass) pairs.

        Parameters
        ----------
        time_s, mass_kg : array_like
            Equal-length 1-D sequences (MKS units).

        Returns
        -------
        numpy.ndarray
            Radiogenic heating [W] for each pair (float64, same length).
        """
        self._check_ptr()
        cdef vector[double] vtime, vmass
        cdef vector[double] vout
        cdef double[::1] mv
        time_c = np.ascontiguousarray(time_s, dtype=np.float64).ravel()
        mass_c = np.ascontiguousarray(mass_kg, dtype=np.float64).ravel()
        mv = time_c; _fill_vector(mv, vtime)
        mv = mass_c; _fill_vector(mv, vmass)
        self._radiogenics_ptr.get().calc_heating_vectorize_all(vtime, vmass, vout)
        return _double_vector_to_ndarray(vout, time_c.shape)

    # ------------------------------------------------------------------------------------------------------------------
    # Config
    # ------------------------------------------------------------------------------------------------------------------
    cpdef dict get_config_dict(self):
        """Return configuration dict with the model name."""
        return PhysicsBase.get_config_dict(self)


# =====================================================================================================================
# OffRadiogenics (alias "none")
# =====================================================================================================================
cdef class OffRadiogenics(RadiogenicsBase):
    """Radiogenics disabled — heating is always zero."""

    def __init__(self):
        cdef c_RadiogenicsConfig config
        cdef c_OffRadiogenics* raw = new c_OffRadiogenics(config)
        self._radiogenics_ptr.reset(<c_RadiogenicsBase*>raw)
        self._ptr = <c_TidalPyBaseClass*>raw


# =====================================================================================================================
# IsotopeRadiogenics
# =====================================================================================================================
cdef class IsotopeRadiogenics(RadiogenicsBase):
    """Radiogenic heating from a set of decaying isotopes.

    Parameters
    ----------
    heat_production_w_kg : array_like
        Per-isotope specific heat production rate [W/kg].
    half_lives_s : array_like
        Per-isotope half life [s].
    mass_fracs : array_like
        Per-isotope mass fraction of the isotope within its element [kg/kg].
    concentrations : array_like
        Per-isotope element concentration in the layer material [kg/kg].
    ref_time_s : float, optional
        Reference time at which the concentrations were measured [s]. Default ``0``.
    names : sequence of str, optional
        Per-isotope labels. Auto-generated as ``isotope_<i>`` when omitted.

    Notes
    -----
    The four numeric arrays are parallel and must have the same length. Heating is
    ``mass * sum_i hpr_i * mass_frac_i * conc_i * exp(ln(0.5) * (t - ref) / t_half_i)``.
    Use :meth:`from_dataset` to build a model from a built-in literature dataset.
    """

    def __cinit__(self, *args, **kwargs):
        self._isotope_ptr = NULL

    def __init__(self, heat_production_w_kg=(), half_lives_s=(),
                 mass_fracs=(), concentrations=(), double ref_time_s=0.0, names=None):
        cdef c_RadiogenicsConfig config
        _build_isotopes(heat_production_w_kg, half_lives_s, mass_fracs,
                        concentrations, names, config.isotopes)
        config.ref_time_s = ref_time_s
        cdef c_IsotopeRadiogenics* raw = new c_IsotopeRadiogenics(config)
        self._radiogenics_ptr.reset(<c_RadiogenicsBase*>raw)
        self._isotope_ptr = raw
        self._ptr         = <c_TidalPyBaseClass*>raw

    def __dealloc__(self):
        self._isotope_ptr = NULL  # base unique_ptr owns the C++ object

    @staticmethod
    def from_dataset(str name):
        """Build an ``IsotopeRadiogenics`` from a built-in literature dataset.

        Parameters
        ----------
        name : str
            A built-in dataset name (see ``available_isotope_datasets()``), e.g.
            ``"modern_day_chondritic"``, ``"llri_and_slri"``, ``"bulk_silicate_earth"``.

        Returns
        -------
        IsotopeRadiogenics
        """
        return make_radiogenics("isotope", {"isotopes": name})

    @property
    def num_isotopes(self) -> int:
        """Number of isotopes in the model."""
        return <int>self._isotope_ptr.get_num_isotopes()

    @property
    def ref_time(self) -> float:
        """Reference time [s]."""
        return self._isotope_ptr.get_ref_time()

    @property
    def isotope_names(self):
        """Per-isotope labels (list of str)."""
        return _isotopes_to_arrays(self._isotope_ptr.get_isotopes())[4]

    @property
    def heat_production(self):
        """Per-isotope specific heat production rate [W/kg]."""
        return _isotopes_to_arrays(self._isotope_ptr.get_isotopes())[0]

    @property
    def half_lives(self):
        """Per-isotope half life [s]."""
        return _isotopes_to_arrays(self._isotope_ptr.get_isotopes())[1]

    @property
    def mass_fracs(self):
        """Per-isotope isotopic mass fraction [kg/kg]."""
        return _isotopes_to_arrays(self._isotope_ptr.get_isotopes())[2]

    @property
    def concentrations(self):
        """Per-isotope element concentration [kg/kg]."""
        return _isotopes_to_arrays(self._isotope_ptr.get_isotopes())[3]

    cpdef dict get_config_dict(self):
        """Return config dict with model name and isotope arrays (as lists)."""
        d = RadiogenicsBase.get_config_dict(self)
        hpr, half, frac, conc, names = _isotopes_to_arrays(self._isotope_ptr.get_isotopes())
        d["heat_production_w_kg"] = list(hpr)
        d["half_lives_s"]         = list(half)
        d["mass_fracs"]           = list(frac)
        d["concentrations"]       = list(conc)
        d["isotope_names"]        = names
        d["ref_time_s"]           = self._isotope_ptr.get_ref_time()
        return d


# =====================================================================================================================
# FixedRadiogenics (alias "constant")
# =====================================================================================================================
cdef class FixedRadiogenics(RadiogenicsBase):
    """Radiogenic heating from a single lumped rate with optional decay.

    Parameters
    ----------
    fixed_heat_production_w_kg : float
        Lumped specific heat production rate [W/kg].
    average_half_life_s : float, optional
        Half life for the lumped rate's exponential decay [s]. Use ``0`` (the
        default) for no decay (a constant heating rate).
    ref_time_s : float, optional
        Reference time at which the rate was measured [s]. Default ``0``.
    """

    def __cinit__(self, *args, **kwargs):
        self._fixed_ptr = NULL

    def __init__(self, double fixed_heat_production_w_kg=0.0,
                 double average_half_life_s=0.0, double ref_time_s=0.0):
        cdef c_RadiogenicsConfig config
        config.fixed_heat_production_w_kg = fixed_heat_production_w_kg
        config.average_half_life_s        = average_half_life_s
        config.ref_time_s                 = ref_time_s
        cdef c_FixedRadiogenics* raw = new c_FixedRadiogenics(config)
        self._radiogenics_ptr.reset(<c_RadiogenicsBase*>raw)
        self._fixed_ptr = raw
        self._ptr       = <c_TidalPyBaseClass*>raw

    def __dealloc__(self):
        self._fixed_ptr = NULL

    @property
    def fixed_heat_production(self) -> float:
        """Lumped specific heat production rate [W/kg]."""
        return self._fixed_ptr.get_fixed_heat_production()

    @property
    def average_half_life(self) -> float:
        """Half life for the lumped rate's decay [s]."""
        return self._fixed_ptr.get_average_half_life()

    @property
    def ref_time(self) -> float:
        """Reference time [s]."""
        return self._fixed_ptr.get_ref_time()

    cpdef dict get_config_dict(self):
        """Return config dict with model name and fixed-rate parameters."""
        d = RadiogenicsBase.get_config_dict(self)
        d["fixed_heat_production_w_kg"] = self._fixed_ptr.get_fixed_heat_production()
        d["average_half_life_s"]        = self._fixed_ptr.get_average_half_life()
        d["ref_time_s"]                 = self._fixed_ptr.get_ref_time()
        return d


# =====================================================================================================================
# Built-in isotope datasets
# =====================================================================================================================
def available_isotope_datasets():
    """Return the list of built-in isotope dataset names.

    These are literature-sourced isotope sets that can be passed to
    ``make_radiogenics("isotope", {"isotopes": <name>})`` or
    ``IsotopeRadiogenics.from_dataset(<name>)``.

    Returns
    -------
    list of str
    """
    return [name.decode("utf-8") for name in c_isotope_dataset_names()]


def isotope_dataset(str name):
    """Return a built-in isotope dataset as an MKS dict.

    Parameters
    ----------
    name : str
        A built-in dataset name (see ``available_isotope_datasets()``).

    Returns
    -------
    dict
        Keys ``heat_production_w_kg``, ``half_lives_s``, ``mass_fracs``,
        ``concentrations`` (lists), ``isotope_names`` (list of str), and
        ``ref_time_s`` (float). All values are MKS.

    Raises
    ------
    ValueError
        If the dataset name is not a built-in (from the C++ catalog).
    """
    cdef c_IsotopeDataset ds = c_get_isotope_dataset(name.encode("utf-8"))
    hpr, half, frac, conc, names = _isotopes_to_arrays(ds.isotopes)
    return {
        "heat_production_w_kg": list(hpr),
        "half_lives_s":         list(half),
        "mass_fracs":           list(frac),
        "concentrations":       list(conc),
        "isotope_names":        names,
        "ref_time_s":           ds.ref_time_s,
    }


# =====================================================================================================================
# Isotope config resolution (non-built-in: explicit arrays, global config, inline dict)
# =====================================================================================================================
def _resolve_isotope_config(dict config):
    """Resolve isotope arrays (in MKS seconds) from a config dict.

    Handles the non-built-in cases (built-in datasets are resolved directly from
    the C++ catalog in ``make_radiogenics``):

    1. Explicit MKS arrays via the keys ``heat_production_w_kg``, ``half_lives_s``,
       ``mass_fracs``, ``concentrations`` (and optional ``ref_time_s``, ``isotope_names``).
    2. A named global-config dataset or inline dict via the ``isotopes`` key. A
       string names a dataset under
       ``TidalPy.config['physics']['radiogenics']['known_isotope_data']``; a dict is
       treated as an inline dataset. These store half lives and reference times in
       Myr, converted to seconds here.

    Returns
    -------
    tuple
        ``(heat_production, half_lives_s, mass_fracs, concentrations, names, ref_time_s)``
        with all numeric lists in MKS units, or ``(None, ...)`` if no isotope data.
    """
    # Explicit MKS arrays take priority.
    if "half_lives_s" in config or "heat_production_w_kg" in config:
        return (
            list(config.get("heat_production_w_kg", ())),
            list(config.get("half_lives_s", ())),
            list(config.get("mass_fracs", ())),
            list(config.get("concentrations", ())),
            config.get("isotope_names", None),
            config.get("ref_time_s", None),
        )

    isotopes = config.get("isotopes", None)
    if isotopes is None:
        return (None, None, None, None, None, None)

    # Resolve a named dataset against the global TidalPy config.
    if isinstance(isotopes, str):
        import TidalPy
        known = TidalPy.config['physics']['radiogenics']['known_isotope_data']
        if isotopes not in known:
            raise ValueError(f"TidalPy: unknown isotope dataset '{isotopes}'.")
        iso_data = known[isotopes]
    elif isinstance(isotopes, dict):
        iso_data = isotopes
    else:
        raise TypeError("TidalPy: 'isotopes' must be a dataset name (str) or an inline dict.")

    # Dataset-level reference time (Myr -> s).
    ref_time_myr = iso_data.get("ref_time", iso_data.get("reference_time", None))
    ref_time_s = None if ref_time_myr is None else ref_time_myr * _SECONDS_PER_MYR

    names = []
    hpr = []
    half_lives_s = []
    mass_fracs = []
    concentrations = []
    for name, entry in iso_data.items():
        if name in ("ref_time", "reference_time"):
            continue
        names.append(name)
        hpr.append(entry["hpr"])
        half_lives_s.append(entry["half_life"] * _SECONDS_PER_MYR)
        mass_fracs.append(entry["iso_mass_fraction"])
        concentrations.append(entry["element_concentration"])

    return (hpr, half_lives_s, mass_fracs, concentrations, names, ref_time_s)


# =====================================================================================================================
# Factory
# =====================================================================================================================
def make_radiogenics(str model_name, dict config=None):
    """Build a radiogenics model from a (case-insensitive) name and config dict.

    This wraps the C++ enum factory: the name (and any alias) is mapped to a
    ``c_RadiogenicsModel`` enum by ``c_radiogenics_model_from_name``, the matching
    concrete model is heap-allocated by ``c_find_radiogenics`` (returning a
    ``unique_ptr``), and the owning pointer is adopted into the correct rich
    Python wrapper.

    Parameters
    ----------
    model_name : str
        Model name or alias. Recognized names: ``off`` (``none``), ``isotope``
        (``isotopes``), ``fixed`` (``constant``).
    config : dict, optional
        Model parameters. For ``isotope``: either explicit MKS arrays
        (``heat_production_w_kg``, ``half_lives_s``, ``mass_fracs``,
        ``concentrations``, ``ref_time_s``) or a named/inline dataset via
        ``isotopes`` (dataset half lives/reference times are in Myr and converted
        to seconds). For ``fixed``: ``fixed_heat_production_w_kg``,
        ``average_half_life_s``, ``ref_time_s``.

    Returns
    -------
    RadiogenicsBase
        A concrete radiogenics model instance.

    Raises
    ------
    ValueError
        If the model name is not recognized.
    """
    if config is None:
        config = {}

    cdef c_RadiogenicsConfig cfg
    cdef c_IsotopeDataset ds

    # Fixed-model scalars (ignored by other models).
    cfg.fixed_heat_production_w_kg = config.get("fixed_heat_production_w_kg", 0.0)
    cfg.average_half_life_s        = config.get("average_half_life_s", 0.0)
    cfg.ref_time_s                 = config.get("ref_time_s", 0.0)

    # Isotope-model data (ignored by other models). A built-in dataset name is
    # resolved straight from the C++ catalog (already MKS); everything else
    # (explicit MKS arrays, global-config dataset, inline dict) goes through the
    # Python resolver (which converts Myr -> s where needed).
    isotopes = config.get("isotopes", None)
    cdef cpp_bool built_in = (
        isinstance(isotopes, str)
        and isotopes.lower() in {name.decode("utf-8") for name in c_isotope_dataset_names()}
    )
    if built_in:
        ds = c_get_isotope_dataset(isotopes.encode("utf-8"))
        cfg.isotopes   = ds.isotopes
        cfg.ref_time_s = ds.ref_time_s
    else:
        hpr, half_lives_s, mass_fracs, concentrations, names, iso_ref = \
            _resolve_isotope_config(config)
        if hpr is not None:
            _build_isotopes(hpr, half_lives_s, mass_fracs, concentrations, names, cfg.isotopes)
            if iso_ref is not None:
                cfg.ref_time_s = iso_ref

    # Map name/alias -> enum (raises ValueError on unknown name via except +),
    # then build the model through the canonical C++ enum factory.
    cdef c_RadiogenicsModel model = c_radiogenics_model_from_name(model_name.encode("utf-8"))
    cdef unique_ptr[c_RadiogenicsBase] ptr = c_find_radiogenics(model, cfg)

    # Adopt the owning unique_ptr into the matching rich Python wrapper.
    cdef OffRadiogenics     o
    cdef IsotopeRadiogenics i
    cdef FixedRadiogenics   f

    if model == c_RadiogenicsModel.Off:
        o = OffRadiogenics.__new__(OffRadiogenics)
        o._radiogenics_ptr = move(ptr)
        o._ptr = <c_TidalPyBaseClass*>o._radiogenics_ptr.get()
        return o
    elif model == c_RadiogenicsModel.Isotope:
        i = IsotopeRadiogenics.__new__(IsotopeRadiogenics)
        i._isotope_ptr     = <c_IsotopeRadiogenics*>ptr.get()
        i._radiogenics_ptr = move(ptr)
        i._ptr             = <c_TidalPyBaseClass*>i._isotope_ptr
        return i
    else:  # c_RadiogenicsModel.Fixed
        f = FixedRadiogenics.__new__(FixedRadiogenics)
        f._fixed_ptr       = <c_FixedRadiogenics*>ptr.get()
        f._radiogenics_ptr = move(ptr)
        f._ptr             = <c_TidalPyBaseClass*>f._fixed_ptr
        return f


# =====================================================================================================================
# Direct heating convenience functions
#
# Each builds a stack-allocated C++ model from its parameters, solves for the
# radiogenic heating, and returns the result (the C++ model is destroyed when the
# function returns). ``time`` and ``mass`` accept Python floats or NumPy arrays
# (broadcast together); model parameters (isotope arrays, fixed rate) are always
# constants. A float result is returned for all-scalar inputs, otherwise a
# float64 ``ndarray``.
# =====================================================================================================================

def off(time, mass):
    """Radiogenic heating for the Off model [W] (always zero). See module notes."""
    cdef c_RadiogenicsConfig cfg
    cdef c_OffRadiogenics model = c_OffRadiogenics(cfg)
    return _solve_heating(<c_RadiogenicsBase*>&model, time, mass)


def isotope(
        time,
        mass,
        heat_production_w_kg=(),
        half_lives_s=(),
        mass_fracs=(),
        concentrations=(),
        double ref_time_s=0.0,
        names=None):
    """Radiogenic heating for the Isotope model [W].

    ``time`` and ``mass`` may be floats or arrays; the isotope arrays are
    constant parameters (all the same length). ``names`` is optional.
    """
    cdef c_RadiogenicsConfig cfg
    _build_isotopes(heat_production_w_kg, half_lives_s, mass_fracs,
                    concentrations, names, cfg.isotopes)
    cfg.ref_time_s = ref_time_s
    cdef c_IsotopeRadiogenics model = c_IsotopeRadiogenics(cfg)
    return _solve_heating(<c_RadiogenicsBase*>&model, time, mass)


def fixed(
        time,
        mass,
        double fixed_heat_production_w_kg=0.0,
        double average_half_life_s=0.0,
        double ref_time_s=0.0):
    """Radiogenic heating for the Fixed model [W] (lumped rate, optional decay)."""
    cdef c_RadiogenicsConfig cfg
    cfg.fixed_heat_production_w_kg = fixed_heat_production_w_kg
    cfg.average_half_life_s        = average_half_life_s
    cfg.ref_time_s                 = ref_time_s
    cdef c_FixedRadiogenics model = c_FixedRadiogenics(cfg)
    return _solve_heating(<c_RadiogenicsBase*>&model, time, mass)
