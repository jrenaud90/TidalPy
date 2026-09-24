# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Cython wrappers for TidalPy's radiogenics models. Each returns total heating [W] from
``calc_heating(time [s], mass [kg])``.

References
----------
- Hussmann and Spohn (2004); Turcotte and Schubert (2001): chondritic isotope data.
- Castillo-Rogez et al. (2007): long- and short-lived radiogenic isotopes.
"""

from libcpp cimport bool as cpp_bool
from libcpp.memory cimport unique_ptr
from libcpp.string cimport string
from libcpp.utility cimport move
from libcpp.vector cimport vector

cimport numpy as cnp

import numpy as np

cnp.import_array()

import TidalPy
from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address, d_SECONDS_PER_MYR
from TidalPy.Utilities_x.classes_x.classes cimport PhysicsBase, c_TidalPyBaseClass
from TidalPy.Utilities_x.classes_x.classes import check_config_keys, factory_defaults

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())


cdef void cy_fill_vector(const double[::1] src, vector[double]& dst) noexcept nogil:
    cdef Py_ssize_t n = src.shape[0]
    cdef Py_ssize_t i
    dst.resize(n)
    for i in range(n):
        dst[i] = src[i]


cdef object cy_double_vector_to_ndarray(vector[double]& src, tuple shape):
    cdef Py_ssize_t n = <Py_ssize_t>src.size()
    cdef Py_ssize_t i
    cdef cnp.ndarray out = np.empty(n, dtype=np.float64)
    cdef double[::1] mv = out
    with nogil:
        for i in range(n):
            mv[i] = src[i]
    return out.reshape(shape)


cdef void cy_build_isotopes(
        object heat_production,
        object half_lives,
        object mass_fracs,
        object concentrations,
        object names,
        vector[c_Isotope]& dst):
    """Build a std::vector[c_Isotope] from parallel sequences; labels default to ``isotope_<i>``."""
    cdef const double[::1] a_hpr  = np.ascontiguousarray(heat_production, dtype=np.float64).ravel()
    cdef const double[::1] a_half = np.ascontiguousarray(half_lives,         dtype=np.float64).ravel()
    cdef const double[::1] a_frac = np.ascontiguousarray(mass_fracs,           dtype=np.float64).ravel()
    cdef const double[::1] a_conc = np.ascontiguousarray(concentrations,       dtype=np.float64).ravel()
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


cdef object cy_isotopes_to_arrays(const vector[c_Isotope]& isotopes):
    """Extract (heat_production, half_lives, mass_fracs, concentrations, names)."""
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
    cdef list names = []
    cdef const c_Isotope* isotope_ptr = NULL
    for i in range(n):
        isotope_ptr = &isotopes[i]
        m_hpr[i]  = isotope_ptr.heat_production
        m_half[i] = isotope_ptr.half_life
        m_frac[i] = isotope_ptr.mass_frac
        m_conc[i] = isotope_ptr.concentration
        names.append(isotope_ptr.name.decode("utf-8"))
    return hpr, half, frac, conc, names


cdef object cy_solve_heating(c_RadiogenicsBase* model, object time, object mass):
    """Dispatch to the most specific C++ vectorized routine for the given input pattern."""
    cdef cpp_bool t_arr = isinstance(time, np.ndarray)
    cdef cpp_bool m_arr = isinstance(mass, np.ndarray)

    cdef vector[double] vtime, vmass
    cdef vector[double] vout
    cdef const double[::1] mv
    cdef const double[::1] mv2
    cdef cnp.ndarray time_arr, mass_arr, t_b, m_b, t_c, m_c
    cdef tuple out_shape
    cdef double scalar_val

    if not (t_arr or m_arr):
        return float(model.calc_heating(<double>time, <double>mass))

    # Time varies; mass constant.
    if t_arr and not m_arr:
        time_arr = np.ascontiguousarray(time, dtype=np.float64)
        out_shape = np.shape(time_arr)
        mv = time_arr.ravel()
        scalar_val = <double>mass  # coercion from Python needs the GIL, so do it before releasing it
        with nogil:
            cy_fill_vector(mv, vtime)
            model.calc_heating_vectorize_time(vtime, scalar_val, vout)
        return cy_double_vector_to_ndarray(vout, out_shape)

    # Mass varies; time constant.
    if m_arr and not t_arr:
        mass_arr = np.ascontiguousarray(mass, dtype=np.float64)
        out_shape = np.shape(mass_arr)
        mv = mass_arr.ravel()
        scalar_val = <double>time
        with nogil:
            cy_fill_vector(mv, vmass)
            model.calc_heating_vectorize_mass(scalar_val, vmass, vout)
        return cy_double_vector_to_ndarray(vout, out_shape)

    # Both vary.
    t_b, m_b = np.broadcast_arrays(
        np.asarray(time, dtype=np.float64),
        np.asarray(mass, dtype=np.float64))
    t_c = np.ascontiguousarray(t_b)
    m_c = np.ascontiguousarray(m_b)
    out_shape = np.shape(t_c)
    mv = t_c.ravel()
    mv2 = m_c.ravel()
    with nogil:
        cy_fill_vector(mv, vtime)
        cy_fill_vector(mv2, vmass)
        model.calc_heating_vectorize_all(vtime, vmass, vout)
    return cy_double_vector_to_ndarray(vout, out_shape)


cdef class RadiogenicsBase(PhysicsBase):
    """Abstract base for radiogenics models; owns the most-derived C++ model object."""

    def __cinit__(self, *args, **kwargs):
        pass  # unique_ptr auto-inits to nullptr; concrete models set it

    def __init__(self, *args, **kwargs):
        raise TypeError(
            "RadiogenicsBase is abstract; instantiate a concrete model "
            "(OffRadiogenics, IsotopeRadiogenics, FixedRadiogenics)."
        )

    def __dealloc__(self):
        self._radiogenics_ptr.reset()
        self._ptr = NULL

    def calc_heating(self, double time, double mass) -> float:
        """Total radiogenic heating [W] for the given time and mass.

        Parameters
        ----------
        time : float
            Elapsed time [s] (shares its zero point with the model reference time).
        mass : float
            Mass of the radiogenic material [kg].

        Returns
        -------
        float
            Radiogenic heating [W].

        Notes
        -----
        Assumes exponential decay from the model's reference time.
        """
        self._check_ptr()
        return self._radiogenics_ptr.get().calc_heating(time, mass)

    def calc_heating_vectorize_time(self, time, double mass):
        """Radiogenic heating over a time sweep at constant mass."""
        self._check_ptr()
        cdef vector[double] vtime
        cdef vector[double] vout
        cdef const double[::1] mv
        cdef cnp.ndarray time_c = np.ascontiguousarray(time, dtype=np.float64).ravel()
        mv = time_c
        with nogil:
            cy_fill_vector(mv, vtime)
            self._radiogenics_ptr.get().calc_heating_vectorize_time(vtime, mass, vout)
        return cy_double_vector_to_ndarray(vout, (time_c.shape[0],))

    def calc_heating_vectorize_mass(self, double time, mass):
        """Radiogenic heating over a mass sweep at constant time."""
        self._check_ptr()
        cdef vector[double] vmass
        cdef vector[double] vout
        cdef const double[::1] mv
        cdef cnp.ndarray mass_c = np.ascontiguousarray(mass, dtype=np.float64).ravel()
        mv = mass_c
        with nogil:
            cy_fill_vector(mv, vmass)
            self._radiogenics_ptr.get().calc_heating_vectorize_mass(time, vmass, vout)
        return cy_double_vector_to_ndarray(vout, (mass_c.shape[0],))

    def calc_heating_vectorize_all(self, time, mass):
        """Radiogenic heating over element-wise (time, mass) pairs of equal length."""
        self._check_ptr()
        cdef vector[double] vtime, vmass
        cdef vector[double] vout
        cdef const double[::1] mv
        cdef const double[::1] mv2
        cdef cnp.ndarray time_c = np.ascontiguousarray(time, dtype=np.float64).ravel()
        cdef cnp.ndarray mass_c = np.ascontiguousarray(mass, dtype=np.float64).ravel()
        mv = time_c
        mv2 = mass_c
        with nogil:
            cy_fill_vector(mv, vtime)
            cy_fill_vector(mv2, vmass)
            self._radiogenics_ptr.get().calc_heating_vectorize_all(vtime, vmass, vout)
        return cy_double_vector_to_ndarray(vout, (time_c.shape[0],))


cdef class OffRadiogenics(RadiogenicsBase):
    """Radiogenics disabled; heating is always zero."""

    def __init__(self):
        cdef c_RadiogenicsConfig config
        cdef unique_ptr[c_RadiogenicsBase] ptr = c_find_radiogenics(c_RadiogenicsModel.Off, config)
        self._radiogenics_ptr = move(ptr)
        self._ptr = <c_TidalPyBaseClass*>self._radiogenics_ptr.get()


cdef class IsotopeRadiogenics(RadiogenicsBase):
    """Radiogenic heating from a set of decaying isotopes.

    Parameters
    ----------
    heat_production : array_like
        Per-isotope specific heat production rate [W/kg].
    half_lives : array_like
        Per-isotope half life [s].
    mass_fracs : array_like
        Per-isotope mass fraction of the isotope within its element [kg/kg].
    concentrations : array_like
        Per-isotope element concentration in the layer material [kg/kg].
    ref_time : float, optional
        Reference time at which the concentrations were measured [s]. Default ``0``.
    names : sequence of str, optional
        Per-isotope labels. Auto-generated as ``isotope_<i>`` when omitted.

    Notes
    -----
    The four numeric arrays are parallel and must share a length. Heating is
    ``mass * sum_i hpr_i * mass_frac_i * conc_i * exp(ln(0.5) * (t - ref) / t_half_i)``.
    """

    def __cinit__(self, *args, **kwargs):
        self._isotope_ptr = NULL

    def __init__(
            self,
            heat_production=(),
            half_lives=(),
            mass_fracs=(),
            concentrations=(),
            double ref_time=0.0,
            names=None):
        cdef c_RadiogenicsConfig config
        cy_build_isotopes(
            heat_production,
            half_lives,
            mass_fracs,
            concentrations,
            names,
            config.isotopes)
        config.ref_time = ref_time
        cdef unique_ptr[c_RadiogenicsBase] ptr = c_find_radiogenics(c_RadiogenicsModel.Isotope, config)
        self._isotope_ptr = <c_IsotopeRadiogenics*>ptr.get()
        self._radiogenics_ptr = move(ptr)
        self._ptr = <c_TidalPyBaseClass*>self._radiogenics_ptr.get()

    def __dealloc__(self):
        self._isotope_ptr = NULL  # RadiogenicsBase._radiogenics_ptr owns the object

    @staticmethod
    def from_dataset(str name):
        """Build an ``IsotopeRadiogenics`` from a built-in dataset (see ``available_isotope_datasets``)."""
        return make_radiogenics("isotope", {"isotopes": name})

    @property
    def num_isotopes(self) -> int:
        """Number of isotopes in the model."""
        self._check_ptr()
        return <int>self._isotope_ptr.get_num_isotopes()

    @property
    def ref_time(self) -> float:
        """Reference time [s]."""
        self._check_ptr()
        return self._isotope_ptr.get_ref_time()

    @property
    def isotope_names(self):
        """Per-isotope labels (list of str)."""
        self._check_ptr()
        return cy_isotopes_to_arrays(self._isotope_ptr.get_isotopes())[4]

    @property
    def heat_production(self):
        """Per-isotope specific heat production rate [W/kg]."""
        self._check_ptr()
        return cy_isotopes_to_arrays(self._isotope_ptr.get_isotopes())[0]

    @property
    def half_lives(self):
        """Per-isotope half life [s]."""
        self._check_ptr()
        return cy_isotopes_to_arrays(self._isotope_ptr.get_isotopes())[1]

    @property
    def mass_fracs(self):
        """Per-isotope isotopic mass fraction [kg/kg]."""
        self._check_ptr()
        return cy_isotopes_to_arrays(self._isotope_ptr.get_isotopes())[2]

    @property
    def concentrations(self):
        """Per-isotope element concentration [kg/kg]."""
        self._check_ptr()
        return cy_isotopes_to_arrays(self._isotope_ptr.get_isotopes())[3]


cdef class FixedRadiogenics(RadiogenicsBase):
    """Radiogenic heating from a single lumped rate with optional decay.

    Parameters
    ----------
    fixed_heat_production : float
        Lumped specific heat production rate [W/kg].
    average_half_life : float, optional
        Half life for the lumped rate's exponential decay [s]; ``0`` (the default) disables the decay.
    ref_time : float, optional
        Reference time at which the rate was measured [s]. Default ``0``.
    """

    def __cinit__(self, *args, **kwargs):
        self._fixed_ptr = NULL

    def __init__(self, double fixed_heat_production=0.0,
                 double average_half_life=0.0, double ref_time=0.0):
        cdef c_RadiogenicsConfig config
        config.fixed_heat_production = fixed_heat_production
        config.average_half_life = average_half_life
        config.ref_time          = ref_time
        cdef unique_ptr[c_RadiogenicsBase] ptr = c_find_radiogenics(c_RadiogenicsModel.Fixed, config)
        self._fixed_ptr = <c_FixedRadiogenics*>ptr.get()
        self._radiogenics_ptr = move(ptr)
        self._ptr = <c_TidalPyBaseClass*>self._radiogenics_ptr.get()

    def __dealloc__(self):
        self._fixed_ptr = NULL

    @property
    def fixed_heat_production(self) -> float:
        """Lumped specific heat production rate [W/kg]."""
        self._check_ptr()
        return self._fixed_ptr.get_fixed_heat_production()

    @property
    def average_half_life(self) -> float:
        """Half life for the lumped rate's decay [s]."""
        self._check_ptr()
        return self._fixed_ptr.get_average_half_life()

    @property
    def ref_time(self) -> float:
        """Reference time [s]."""
        self._check_ptr()
        return self._fixed_ptr.get_ref_time()


def available_isotope_datasets():
    """Names of the built-in, literature-sourced isotope datasets."""
    return [name.decode("utf-8") for name in c_isotope_dataset_names()]


def isotope_dataset(str name):
    """Return a built-in isotope dataset as an MKS dict, keyed as ``make_radiogenics`` accepts.

    Parameters
    ----------
    name : str
        A built-in dataset name (see ``available_isotope_datasets()``).

    Returns
    -------
    dict

    Raises
    ------
    ValueError
        Unknown dataset name.
    """
    cdef c_IsotopeDataset ds = c_get_isotope_dataset(name.encode("utf-8"))
    hpr, half, frac, conc, names = cy_isotopes_to_arrays(ds.isotopes)
    return {
        "heat_production_w_kg": list(hpr),
        "half_lives_s":         list(half),
        "mass_fracs":           list(frac),
        "concentrations":       list(conc),
        "isotope_names":        names,
        "ref_time_s":           ds.ref_time,
    }


def _resolve_isotope_config(dict config):
    """Resolve isotope arrays in MKS from a config dict.

    Covers what the C++ catalog does not: explicit MKS arrays, or an ``isotopes`` key naming a dataset
    under ``TidalPy.config['physics']['radiogenics']['known_isotope_data']`` or holding an inline dict.
    Those datasets store half lives and reference times in Myr, converted to seconds here.

    Returns
    -------
    tuple
        ``(heat_production, half_lives, mass_fracs, concentrations, names, ref_time)``, all ``None`` when
        the config carries no isotope data.
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

    # `isotopes` is a dataset name or an inline table, so it stays `object` until the isinstance checks
    # below say which. The reference time is a float or absent, hence `object` rather than `double`.
    cdef object isotopes = config.get("isotopes", None)
    cdef dict known
    cdef dict iso_data
    cdef object ref_time_myr
    cdef object ref_time
    cdef str name
    # `object`, not `dict`: the table also carries a scalar `ref_time`/`reference_time`, and the loop below
    # unpacks every item before the `continue` that skips those keys, so `dict` would raise first.
    cdef object entry
    if isotopes is None:
        return (None, None, None, None, None, None)

    if isinstance(isotopes, str):
        known = TidalPy.config['physics']['radiogenics']['known_isotope_data']
        if isotopes not in known:
            raise ValueError(f"TidalPy: unknown isotope dataset '{isotopes}'.")
        iso_data = known[isotopes]
    elif isinstance(isotopes, dict):
        iso_data = isotopes
    else:
        raise TypeError("TidalPy: 'isotopes' must be a dataset name (str) or an inline dict.")

    # Myr -> s.
    ref_time_myr = iso_data.get("ref_time", iso_data.get("reference_time", None))
    ref_time = None if ref_time_myr is None else ref_time_myr * d_SECONDS_PER_MYR

    cdef list names = []
    cdef list hpr = []
    cdef list half_lives = []
    cdef list mass_fracs = []
    cdef list concentrations = []
    for name, entry in iso_data.items():
        if name in ("ref_time", "reference_time"):
            continue
        names.append(name)
        hpr.append(entry["hpr"])
        half_lives.append(entry["half_life"] * d_SECONDS_PER_MYR)
        mass_fracs.append(entry["iso_mass_fraction"])
        concentrations.append(entry["element_concentration"])

    return (hpr, half_lives, mass_fracs, concentrations, names, ref_time)


# Every config key any radiogenics model reads; make_radiogenics rejects anything else.
RADIOGENICS_CONFIG_KEYS = frozenset({
    "fixed_heat_production_w_kg", "average_half_life_s", "ref_time_s", "isotopes",
    "heat_production_w_kg", "half_lives_s", "mass_fracs", "concentrations", "isotope_names"})


def _same_model(str table_name, str model_name) -> bool:
    """Whether two names (aliases included) resolve to the same model."""
    return c_radiogenics_model_from_name(table_name.lower().encode("utf-8")) == c_radiogenics_model_from_name(model_name.lower().encode("utf-8"))


def make_radiogenics(str model_name, dict config=None):
    """Build a radiogenics model from a (case-insensitive) name and config dict.

    Parameters
    ----------
    model_name : str
        Model name or alias: ``off`` (``none``), ``isotope`` (``isotopes``), ``fixed`` (``constant``).
    config : dict, optional
        Model parameters under the unit-suffixed keys ``get_config_dict()`` emits. For ``isotope``: either
        explicit MKS arrays, or a named or inline dataset under ``isotopes`` whose half lives and reference
        times are in Myr and converted here. For ``fixed``: ``fixed_heat_production_w_kg``,
        ``average_half_life_s``, ``ref_time_s``. Each model ignores the other's keys, so a table merged
        family-wide builds the model it names from that model's keys alone.

    Returns
    -------
    RadiogenicsBase

    Raises
    ------
    ValueError
        Unknown model name, or a config key that no radiogenics model reads.
    """
    if config is None:
        # Fall back to the same defaults the world-attached path uses.
        config = factory_defaults("radiogenics", RADIOGENICS_CONFIG_KEYS, model_name, _same_model)
    check_config_keys(config, RADIOGENICS_CONFIG_KEYS, "radiogenics")
    if config is None:
        config = {}

    cdef c_RadiogenicsConfig cfg
    cdef c_IsotopeDataset ds

    cdef c_RadiogenicsModel model = c_radiogenics_model_from_name(model_name.encode("utf-8"))

    # The default-constructed config carries the C++ defaults, so only override what the caller gave.
    if "fixed_heat_production_w_kg" in config:
        cfg.fixed_heat_production = config["fixed_heat_production_w_kg"]
    if "average_half_life_s" in config:
        cfg.average_half_life = config["average_half_life_s"]
    if "ref_time_s" in config:
        cfg.ref_time = config["ref_time_s"]

    # Read only for the isotope model: a family-wide merged config can carry a dataset next to a fixed
    # model, and the dataset's reference time must not become the fixed rate's. A built-in name resolves
    # straight from the C++ catalog (already MKS); anything else goes through the Python resolver.
    cdef object isotopes = config.get("isotopes", None)
    cdef cpp_bool built_in = (
        isinstance(isotopes, str)
        and isotopes.lower() in {name.decode("utf-8") for name in c_isotope_dataset_names()}
    )
    # Explicit MKS arrays win over a named dataset (the resolver's rule), so a dataset name merged in from the
    # material defaults never replaces the isotopes a layer lists itself.
    cdef cpp_bool explicit_arrays = ("half_lives_s" in config) or ("heat_production_w_kg" in config)
    if model == c_RadiogenicsModel.Isotope:
        if built_in and not explicit_arrays:
            ds = c_get_isotope_dataset(isotopes.encode("utf-8"))
            cfg.isotopes = ds.isotopes
            # A given reference time says when the dataset's abundances apply; otherwise the dataset's own.
            cfg.ref_time = config["ref_time_s"] if "ref_time_s" in config else ds.ref_time
        else:
            hpr, half_lives, mass_fracs, concentrations, names, iso_ref = \
                _resolve_isotope_config(config)
            if hpr is not None:
                cy_build_isotopes(hpr, half_lives, mass_fracs, concentrations, names, cfg.isotopes)
                if iso_ref is not None:
                    cfg.ref_time = iso_ref

    cdef unique_ptr[c_RadiogenicsBase] ptr = c_find_radiogenics(model, cfg)

    # Adopt the owning unique_ptr into the matching Python wrapper.
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


# Convenience functions. Each builds a stack-allocated C++ model that dies with the call. ``time`` and
# ``mass`` accept floats or ndarrays broadcast together; the model parameters stay scalar.

def off(time, mass):
    """Radiogenic heating for the Off model [W]; always zero."""
    cdef c_RadiogenicsConfig cfg
    cdef c_OffRadiogenics model = c_OffRadiogenics(cfg)
    return cy_solve_heating(<c_RadiogenicsBase*>&model, time, mass)


def isotope(
        time,
        mass,
        heat_production=(),
        half_lives=(),
        mass_fracs=(),
        concentrations=(),
        double ref_time=0.0,
        names=None):
    """Radiogenic heating for the Isotope model [W]. The isotope arrays must all share a length."""
    cdef c_RadiogenicsConfig cfg
    cy_build_isotopes(
        heat_production,
        half_lives,
        mass_fracs,
        concentrations,
        names,
        cfg.isotopes)
    cfg.ref_time = ref_time
    cdef c_IsotopeRadiogenics model = c_IsotopeRadiogenics(cfg)
    return cy_solve_heating(<c_RadiogenicsBase*>&model, time, mass)


def fixed(
        time,
        mass,
        double fixed_heat_production=0.0,
        double average_half_life=0.0,
        double ref_time=0.0):
    """Radiogenic heating for the Fixed model [W] (lumped rate, optional decay)."""
    cdef c_RadiogenicsConfig cfg
    cfg.fixed_heat_production = fixed_heat_production
    cfg.average_half_life = average_half_life
    cfg.ref_time          = ref_time
    cdef c_FixedRadiogenics model = c_FixedRadiogenics(cfg)
    return cy_solve_heating(<c_RadiogenicsBase*>&model, time, mass)
