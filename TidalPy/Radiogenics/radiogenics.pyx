# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Cython wrappers for TidalPy's radiogenics models. Each returns total heating [W] from
``calc_heating(time [s], mass [kg])``.

Each model's parameters, defaults, bounds, and descriptions come from its C++ parameter table. The isotope model also
takes its isotopes from a dataset (``isotopes``: a built-in name, a dataset of ``[radiogenics.known_isotope_data]``,
or an inline table) and labels them (``isotope_names``).

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

# The shared vector helpers build their result arrays through the NumPy C API.
cnp.import_array()

import TidalPy
from TidalPy.Utilities.logging.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address, d_SECONDS_PER_MYR
from TidalPy.Utilities.arrays.vectors cimport cy_broadcast_inputs, cy_vector_to_ndarray
from TidalPy.Utilities.classes.classes cimport (
    PhysicsBase,
    c_ParamMap,
    c_share_physics,
    cy_collect_parameters,
    cy_param_map,
)
from TidalPy.Utilities.classes.classes import canonical_parameter_keys, factory_defaults
from TidalPy.Utilities.classes.families import ModelFamily

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())

# The isotope model's per-isotope tables, by config key; a model given any of them takes its isotopes from them rather
# than from a dataset.
_ISOTOPE_TABLE_KEYS = ("heat_production_w_kg", "half_lives_s", "mass_fracs", "concentrations")


cdef object cy_solve_heating(c_RadiogenicsBase* model, object time, object mass, cpp_bool flatten):
    """Heating for float or ndarray ``time`` and ``mass`` broadcast together (cy_broadcast_inputs); a float when both
    are floats and ``flatten`` is off."""
    cdef vector[vector[double]] inputs
    cdef vector[double] heating
    cdef object shape = cy_broadcast_inputs((time, mass), inputs, flatten)
    if shape is None:
        return float(model.calc_heating(<double>time, <double>mass))
    with nogil:
        model.calc_heating_vectorize(inputs[0], inputs[1], heating)
    return cy_vector_to_ndarray(heating, shape)


cdef dict cy_dataset_parameters(const c_IsotopeDataset& dataset):
    """A built-in dataset as the isotope model's parameters (MKS), with its labels under ``isotope_names``."""
    cdef dict parameters = {key: [] for key in _ISOTOPE_TABLE_KEYS}
    cdef list names = []
    cdef size_t isotope_i
    for isotope_i in range(dataset.isotopes.size()):
        parameters["heat_production_w_kg"].append(dataset.isotopes[isotope_i].heat_production)
        parameters["half_lives_s"].append(dataset.isotopes[isotope_i].half_life)
        parameters["mass_fracs"].append(dataset.isotopes[isotope_i].mass_frac)
        parameters["concentrations"].append(dataset.isotopes[isotope_i].concentration)
        names.append(dataset.isotopes[isotope_i].name.decode("utf-8"))
    parameters["ref_time_s"] = dataset.ref_time
    parameters["isotope_names"] = names
    return parameters


cdef class RadiogenicsBase(PhysicsBase):
    """Base for radiogenics models, which give a layer's radiogenic heating (``calc_heating``).

    Instantiate a concrete model (``OffRadiogenics``, ``IsotopeRadiogenics``, ``FixedRadiogenics``) with its
    parameters positionally (in the order ``get_parameter_info()`` lists them) or as keywords (argument names or config
    keys), or build one by name with ``make_radiogenics``. A layer shares the model it is given
    (``Layer.radiogenics``).
    """

    # The canonical name of the model a concrete subclass builds; None on this abstract base.
    MODEL_NAME = None

    def __init__(self, *args, dict config=None, **parameters):
        cdef object model_name = type(self).MODEL_NAME
        if model_name is None:
            raise TypeError(
                "RadiogenicsBase is abstract; instantiate a concrete model (OffRadiogenics, IsotopeRadiogenics, "
                "FixedRadiogenics) or call make_radiogenics.")
        cdef c_ParamMap param_map = cy_param_map(cy_collect_parameters(type(self), args, config, parameters))
        cdef unique_ptr[c_RadiogenicsBase] model = c_find_radiogenics((<str>model_name).encode("utf-8"), param_map)
        self._set_model(c_share_physics[c_RadiogenicsBase](move(model)))

    cdef c_RadiogenicsBase* _radiogenics(self) except NULL:
        self._check_ptr()
        return <c_RadiogenicsBase*>self._ptr

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
        return self._radiogenics().calc_heating(time, mass)

    def calc_heating_vectorize_time(self, time, double mass):
        """Radiogenic heating over a time sweep at constant mass."""
        return cy_solve_heating(self._radiogenics(), time, mass, True)

    def calc_heating_vectorize_mass(self, double time, mass):
        """Radiogenic heating over a mass sweep at constant time."""
        return cy_solve_heating(self._radiogenics(), time, mass, True)

    def calc_heating_vectorize_all(self, time, mass):
        """Radiogenic heating over element-wise (time, mass) pairs of equal length."""
        return cy_solve_heating(self._radiogenics(), time, mass, True)


cdef class OffRadiogenics(RadiogenicsBase):
    """Radiogenics disabled (alias ``"none"``); heating is always zero."""
    MODEL_NAME = "off"


cdef class IsotopeRadiogenics(RadiogenicsBase):
    """Radiogenic heating from a set of decaying isotopes (alias ``"isotopes"``):
    ``mass * sum_i hpr_i * mass_frac_i * conc_i * exp(ln(0.5) * (t - ref_time) / t_half_i)``.

    The isotopes come from the four parallel tables (``heat_production``, ``half_lives``, ``mass_fracs``,
    ``concentrations``, one value per isotope), or, when none of them is given, from a dataset: ``isotopes`` names a
    built-in one (``available_isotope_datasets()``) or one of ``[radiogenics.known_isotope_data]``, or holds an
    inline table in the same form; with neither, the ``[radiogenics] isotopes`` dataset of the TidalPy configuration.
    A dataset brings its labels and reference time; a given ``ref_time`` says when its abundances apply instead.
    ``isotope_names`` labels the isotopes, one per isotope; labels take no part in the heating.
    """
    MODEL_NAME = "isotope"
    # Keys the isotope model reads beyond its parameter table.
    EXTRA_CONFIG_KEYS = ("isotopes", "isotope_names")

    def __init__(self, *args, dict config=None, **parameters):
        cdef dict merged = cy_collect_parameters(type(self), args, config, parameters)
        cdef object names = merged.pop("isotope_names", None)
        cdef object dataset = merged.pop("isotopes", None)
        cdef dict from_dataset
        if not any(merged.get(key) is not None for key in _ISOTOPE_TABLE_KEYS):
            if dataset is None:
                dataset = factory_defaults("radiogenics", ("isotopes",)).get("isotopes")
            if dataset is not None:
                from_dataset = isotope_dataset_parameters(dataset)
                if names is None:
                    names = from_dataset["isotope_names"]
                del from_dataset["isotope_names"]
                merged = {**from_dataset, **{key: value for key, value in merged.items() if value is not None}}
        if isinstance(names, str):
            raise TypeError("TidalPy: 'isotope_names' takes one name per isotope, not a single string.")
        cdef vector[string] c_names
        if names is not None:
            for name in names:
                c_names.push_back(str(name).encode("utf-8"))
        cdef unique_ptr[c_RadiogenicsBase] model = c_make_isotope_radiogenics(cy_param_map(merged), c_names)
        self._set_model(c_share_physics[c_RadiogenicsBase](move(model)))

    @classmethod
    def _default_model(cls):
        """The model with no isotopes, which reads no configuration (for its parameter descriptions)."""
        cdef IsotopeRadiogenics model = cls.__new__(cls)
        cdef c_ParamMap no_parameters
        cdef unique_ptr[c_RadiogenicsBase] built = c_find_radiogenics(b"isotope", no_parameters)
        model._set_model(c_share_physics[c_RadiogenicsBase](move(built)))
        return model

    def with_parameters(self, **changes):
        """A new model with the given parameters (argument names or config keys) and ``isotope_names`` changed.

        Raises
        ------
        ValueError
            An unknown parameter, a value outside its bounds, tables of different lengths, or labels that do not
            match the isotopes.
        """
        if "isotope_names" not in changes:
            return RadiogenicsBase.with_parameters(self, **changes)
        cdef dict current = {key: value for key, value in self.get_config_dict().items() if key != "model"}
        return type(self)(config={**current, **canonical_parameter_keys(type(self), changes)})

    @property
    def num_isotopes(self) -> int:
        """Number of isotopes in the model."""
        return <int>(<c_IsotopeRadiogenics*>self._radiogenics()).get_num_isotopes()

    @property
    def isotope_names(self) -> list:
        """The isotopes' labels, in isotope order; empty for unlabeled isotopes."""
        cdef vector[string] names = (<c_IsotopeRadiogenics*>self._radiogenics()).get_isotope_names()
        return [name.decode("utf-8") for name in names]


cdef class FixedRadiogenics(RadiogenicsBase):
    """Radiogenic heating from one lumped rate with optional decay (alias ``"constant"``):
    ``mass * rate * exp(ln(0.5) * (t - ref_time) / average_half_life)``, with no decay for an average half life of 0
    or below."""
    MODEL_NAME = "fixed"


def available_isotope_datasets():
    """Names of the built-in, literature-sourced isotope datasets."""
    return [name.decode("utf-8") for name in c_isotope_dataset_names()]


def isotope_dataset(str name):
    """A built-in isotope dataset as the isotope model's parameters (MKS), keyed as ``make_radiogenics`` takes them.

    Parameters
    ----------
    name : str
        A built-in dataset name, case-insensitive (see ``available_isotope_datasets()``).

    Returns
    -------
    dict
        The four per-isotope tables, ``ref_time_s``, and ``isotope_names``.

    Raises
    ------
    ValueError
        Unknown dataset name.
    """
    return cy_dataset_parameters(c_get_isotope_dataset(name.encode("utf-8")))


def isotope_dataset_parameters(object dataset) -> dict:
    """An isotope dataset as the isotope model's parameters (MKS), with its labels under ``isotope_names``.

    Parameters
    ----------
    dataset : str or dict
        A built-in dataset name (case-insensitive), a dataset of ``[radiogenics.known_isotope_data]`` in the TidalPy
        configuration, or an inline table in that form: a ``ref_time`` (or ``reference_time``) in Myr and one table
        per isotope holding ``hpr`` [W/kg], ``half_life`` [Myr], ``iso_mass_fraction``, and
        ``element_concentration``. A built-in name wins over a configured dataset of the same name.

    Raises
    ------
    ValueError
        Unknown dataset name.
    TypeError
        ``dataset`` is neither a name nor a table.
    """
    cdef dict table
    cdef dict parameters
    cdef object ref_time_myr
    if isinstance(dataset, str):
        if dataset.lower() in available_isotope_datasets():
            return isotope_dataset(dataset)
        known = ((TidalPy.config or {}).get("radiogenics", {}) or {}).get("known_isotope_data", {}) or {}
        if dataset not in known:
            raise ValueError(
                f"TidalPy: unknown isotope dataset '{dataset}'. Built in: {', '.join(available_isotope_datasets())}; "
                f"configured: {', '.join(known) if known else 'none'}.")
        table = known[dataset]
    elif isinstance(dataset, dict):
        table = dataset
    else:
        raise TypeError("TidalPy: 'isotopes' is a dataset name or an inline dataset table.")

    parameters = {key: [] for key in _ISOTOPE_TABLE_KEYS}
    parameters["isotope_names"] = []
    for name, entry in table.items():
        if name in ("ref_time", "reference_time"):
            continue
        parameters["isotope_names"].append(name)
        parameters["heat_production_w_kg"].append(entry["hpr"])
        parameters["half_lives_s"].append(entry["half_life"] * d_SECONDS_PER_MYR)
        parameters["mass_fracs"].append(entry["iso_mass_fraction"])
        parameters["concentrations"].append(entry["element_concentration"])
    ref_time_myr = table.get("ref_time", table.get("reference_time"))
    if ref_time_myr is not None:
        parameters["ref_time_s"] = ref_time_myr * d_SECONDS_PER_MYR
    return parameters


def _canonical_name(str model_name) -> str:
    return c_radiogenics_canonical_name(model_name.encode("utf-8")).decode("utf-8")


_FAMILY = ModelFamily("radiogenics", (OffRadiogenics, IsotopeRadiogenics, FixedRadiogenics), _canonical_name)

# Every config key any radiogenics model reads.
RADIOGENICS_CONFIG_KEYS = _FAMILY.config_keys


def radiogenics_model_names() -> tuple:
    """The canonical names of the radiogenics models."""
    return _FAMILY.model_names()


def radiogenics_config_keys(str model_name) -> frozenset:
    """The config keys a radiogenics model reads, by any of its names.

    Raises
    ------
    ValueError
        Unknown model name.
    """
    return _FAMILY.config_keys_of(model_name)


def make_radiogenics(str model_name, dict config=None):
    """Build a radiogenics model from a (case-insensitive) name and config dict.

    Parameters
    ----------
    model_name : str
        Model name or alias: ``off`` (``none``), ``isotope`` (``isotopes``), ``fixed`` (``constant``).
    config : dict, optional
        Model parameters by config key (see each model's ``get_parameter_info()``); absent keys (all of them for
        ``None``) take the model's defaults. The isotope model also reads ``isotopes`` and ``isotope_names`` (see
        ``IsotopeRadiogenics``), and one given no isotopes takes the ``[radiogenics] isotopes`` dataset of the TidalPy
        configuration.

    Returns
    -------
    RadiogenicsBase

    Raises
    ------
    ValueError
        Unknown model name, or a parameter the model does not read; each message names the closest accepted one.
    """
    return _FAMILY.make(model_name, config)


# Convenience functions. Each builds the model for the one call. ``time`` and ``mass`` accept floats or ndarrays
# broadcast together; the model parameters stay scalar or, for the isotope tables, one value per isotope.

cdef object cy_direct_heating(str model_name, dict parameters, object time, object mass):
    """Heating from a model built for the one call (by name and parameters), over ``time`` and ``mass``."""
    cdef unique_ptr[c_RadiogenicsBase] model = c_find_radiogenics(model_name.encode("utf-8"), cy_param_map(parameters))
    return cy_solve_heating(model.get(), time, mass, False)


def off(time, mass):
    """Radiogenic heating for the Off model [W]; always zero."""
    return cy_direct_heating("off", {}, time, mass)


def isotope(
        time,
        mass,
        heat_production=(),
        half_lives=(),
        mass_fracs=(),
        concentrations=(),
        double ref_time=0.0):
    """Radiogenic heating for the Isotope model [W]. The four isotope tables hold one value per isotope."""
    return cy_direct_heating(
        "isotope",
        {"heat_production": heat_production, "half_lives": half_lives, "mass_fracs": mass_fracs,
         "concentrations": concentrations, "ref_time": ref_time},
        time,
        mass)


def fixed(
        time,
        mass,
        double fixed_heat_production=0.0,
        double average_half_life=0.0,
        double ref_time=0.0):
    """Radiogenic heating for the Fixed model [W] (lumped rate, optional decay)."""
    return cy_direct_heating(
        "fixed",
        {"fixed_heat_production": fixed_heat_production, "average_half_life": average_half_life,
         "ref_time": ref_time},
        time,
        mass)
