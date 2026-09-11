# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""Native input builders for the `_x` radial solver.

`build_rs_input_homogeneous_layers` and `build_rs_input_from_data` assemble the array-based inputs
that :func:`TidalPy.RadialSolver_x.radial_solver` expects from a per-layer description of the planet.

Grid repairs made by the from-data builder are logged as warnings through the TidalPy logger.
"""

from libc.string cimport memcpy
from libcpp cimport bool as cpp_bool
from libcpp.complex cimport complex as cpp_complex
from libcpp.vector cimport vector

from collections import namedtuple

import numpy as np
cimport numpy as cnp
cnp.import_array()

from TidalPy.Utilities_x.logging_x.logger cimport (
    set_tidalpy_logger_ptr_void,
    get_tidalpy_logger_address,
)
from TidalPy.constants cimport set_tidalpy_config_ptr, get_shared_config_address
from TidalPy.rheology_x.rheology cimport RheologyBase, c_RheologyBase
from TidalPy.rheology_x.rheology import make_rheology

# Wire this DLL's shared pointers to the process-wide TidalPy singletons.
set_tidalpy_logger_ptr_void(get_tidalpy_logger_address())
set_tidalpy_config_ptr(get_shared_config_address())


PlanetBuildData = namedtuple(
    "PlanetBuildData",
    (
        "radius_array",
        "density_array",
        "complex_bulk_modulus_array",
        "complex_shear_modulus_array",
        "frequency",
        "planet_bulk_density",
        "layer_types",
        "is_static_bylayer",
        "is_incompressible_bylayer",
        "upper_radius_bylayer_array",
    ))
PlanetBuildData.__doc__ = (
    "Radial-solver inputs assembled by the `_x` builders. The fields are, in order, the positional "
    "arguments of `TidalPy.RadialSolver_x.radial_solver`, so ``radial_solver(*build_data, **kwargs)`` "
    "runs the solve.")


# =====================================================================================================================
# Internal conversion helpers
# =====================================================================================================================

cdef tuple _as_layer_tuple(object values, size_t num_layers, str name):
    """Return `values` as a tuple with one entry per layer (ValueError otherwise)."""
    if values is None:
        raise ValueError(f"`{name}` must be provided with one entry per layer.")
    cdef tuple out = tuple(values)
    if <size_t>len(out) != num_layers:
        raise ValueError(
            f"`{name}` must have one entry per layer ({num_layers}), found {len(out)}.")
    return out


cdef void _fill_double_vector(object values, str name, vector[double]& out) except *:
    """Copy a sequence of real numbers into a C++ vector."""
    cdef size_t i
    cdef size_t n = <size_t>len(values)
    out.clear()
    out.reserve(n)
    try:
        for i in range(n):
            out.push_back(<double>values[i])
    except (TypeError, ValueError) as exc:
        raise ValueError(f"`{name}` must contain real numbers; entry {i} is {values[i]!r}.") from exc


cdef void _fill_size_vector(object values, str name, vector[size_t]& out) except *:
    """Copy a sequence of non-negative integers into a C++ vector."""
    cdef size_t i
    cdef size_t n = <size_t>len(values)
    cdef object value
    out.clear()
    out.reserve(n)
    for i in range(n):
        value = values[i]
        if not isinstance(value, (int, np.integer)) or value < 0:
            raise ValueError(f"`{name}` must contain non-negative integers; entry {i} is {value!r}.")
        out.push_back(<size_t>value)


cdef object _as_f64_1d(object values, str name):
    """Return `values` as a C-contiguous 1-D float64 array."""
    arr = np.ascontiguousarray(values, dtype=np.float64)
    if arr.ndim != 1:
        raise ValueError(f"`{name}` must be one dimensional, found {arr.ndim} dimensions.")
    return arr


cdef void _fill_vector_from_array(object values, str name, vector[double]& out) except *:
    """Copy a 1-D float64 array into a C++ vector."""
    cdef cnp.ndarray[cnp.float64_t, ndim=1] arr = _as_f64_1d(values, name)
    cdef size_t n = <size_t>arr.shape[0]
    out.resize(n)
    if n > 0:
        memcpy(out.data(), cnp.PyArray_DATA(arr), n * sizeof(double))


cdef cnp.ndarray _vector_to_f64_array(const vector[double]& vec):
    """Copy a C++ double vector into a new float64 numpy array."""
    cdef cnp.npy_intp n = <cnp.npy_intp>vec.size()
    cdef cnp.ndarray[cnp.float64_t, ndim=1] arr = np.empty(n, dtype=np.float64, order="C")
    if n > 0:
        memcpy(cnp.PyArray_DATA(arr), vec.data(), <size_t>n * sizeof(double))
    return arr


cdef cnp.ndarray _vector_to_c128_array(const vector[cpp_complex[double]]& vec):
    """Copy a C++ complex vector into a new complex128 numpy array."""
    cdef cnp.npy_intp n = <cnp.npy_intp>vec.size()
    cdef cnp.ndarray[cnp.complex128_t, ndim=1] arr = np.empty(n, dtype=np.complex128, order="C")
    if n > 0:
        memcpy(cnp.PyArray_DATA(arr), vec.data(), <size_t>n * sizeof(cpp_complex[double]))
    return arr


cdef object _coerce_rheology(object model, str argument_name, object position):
    """Return a `rheology_x` model instance for `model` (an instance or a model name)."""
    if isinstance(model, RheologyBase):
        return model
    if isinstance(model, str):
        return make_rheology(model)
    where = f"`{argument_name}`" if position is None else f"`{argument_name}` entry {position}"
    hint = ""
    if type(model).__module__.startswith("TidalPy.rheology."):
        hint = (" Classic `TidalPy.rheology` models are not accepted by the `_x` builders; use the "
                "matching model from `TidalPy.rheology_x`.")
    raise TypeError(
        f"{where} must be a `TidalPy.rheology_x` model instance or a model name, "
        f"found {type(model).__name__}.{hint}")


cdef list cy_resolve_rheology_bylayer(
        object models, size_t num_layers, str argument_name, vector[const c_RheologyBase*]& out_ptrs):
    """Resolve a rheology argument into one C++ model pointer per layer.

    `models` may be a single model instance or name (applied to every layer) or a sequence with one
    instance or name per layer. The returned list holds the Python objects that own the pointers.
    """
    cdef list keep_alive = []
    cdef size_t layer_i
    cdef RheologyBase wrapper
    cdef const c_RheologyBase* model_ptr

    out_ptrs.clear()
    out_ptrs.reserve(num_layers)
    if isinstance(models, (RheologyBase, str)):
        wrapper = <RheologyBase>_coerce_rheology(models, argument_name, None)
        keep_alive.append(wrapper)
        model_ptr = wrapper._rheology_ptr.get()
        if model_ptr == NULL:
            raise ValueError(f"`{argument_name}` model is not initialized.")
        for layer_i in range(num_layers):
            out_ptrs.push_back(model_ptr)
        return keep_alive

    if models is None:
        raise ValueError(
            f"`{argument_name}` must be a `TidalPy.rheology_x` model (applied to every layer) or a "
            f"sequence with one model per layer.")
    try:
        models_seq = tuple(models)
    except TypeError:
        _coerce_rheology(models, argument_name, None)  # Raises the informative TypeError.
    if <size_t>len(models_seq) != num_layers:
        raise ValueError(
            f"`{argument_name}` must have one rheology model per layer ({num_layers}), "
            f"found {len(models_seq)}. Pass a single model to apply it to every layer.")
    for layer_i in range(num_layers):
        wrapper = <RheologyBase>_coerce_rheology(models_seq[layer_i], argument_name, layer_i)
        keep_alive.append(wrapper)
        model_ptr = wrapper._rheology_ptr.get()
        if model_ptr == NULL:
            raise ValueError(f"`{argument_name}` entry {layer_i} is not initialized.")
        out_ptrs.push_back(model_ptr)
    return keep_alive


cdef object _build_outputs(
        const c_RadialSolverInputs& inputs, tuple layer_types, tuple is_static_bylayer,
        tuple is_incompressible_bylayer):
    """Convert the C++ result into a `PlanetBuildData` namedtuple of numpy arrays."""
    return PlanetBuildData(
        _vector_to_f64_array(inputs.radius_m),
        _vector_to_f64_array(inputs.density_kg_m3),
        _vector_to_c128_array(inputs.complex_bulk_modulus_pa),
        _vector_to_c128_array(inputs.complex_shear_modulus_pa),
        inputs.forcing_frequency_rad_s,
        inputs.planet_bulk_density_kg_m3,
        layer_types,
        is_static_bylayer,
        is_incompressible_bylayer,
        _vector_to_f64_array(inputs.upper_radius_bylayer_m),
    )


# =====================================================================================================================
# Public builders
# =====================================================================================================================

def build_rs_input_homogeneous_layers(
        double planet_radius,
        double forcing_frequency,
        density_tuple,
        static_bulk_modulus_tuple,
        static_shear_modulus_tuple,
        bulk_viscosity_tuple,
        shear_viscosity_tuple,
        layer_type_tuple,
        layer_is_static_tuple,
        layer_is_incompressible_tuple,
        shear_rheology_model_tuple,
        bulk_rheology_model_tuple,
        radius_fraction_tuple=None,
        thickness_fraction_tuple=None,
        volume_fraction_tuple=None,
        slices_tuple=None,
        size_t slice_per_layer=10,
        cpp_bool perform_checks=True):
    """Build radial-solver inputs for a planet whose layers each have constant properties.

    The radial grid of each layer runs from the layer base to its top (both inclusive) with
    `slices_tuple[i]` (or `slice_per_layer`) evenly spaced slices, so interface radii appear twice as
    the solver requires. The complex moduli are evaluated at `forcing_frequency` with the supplied
    ``rheology_x`` models.

    Parameters
    ----------
    planet_radius : float
        Planet radius [m].
    forcing_frequency : float
        Tidal forcing frequency [rad s-1].
    density_tuple, static_bulk_modulus_tuple, static_shear_modulus_tuple : sequence of float
        Per-layer density [kg m-3] and unrelaxed bulk and shear moduli [Pa].
    bulk_viscosity_tuple, shear_viscosity_tuple : sequence of float
        Per-layer bulk and shear viscosities [Pa s].
    layer_type_tuple : sequence of str
        Per-layer type ("solid" or "liquid"); passed through to the output.
    layer_is_static_tuple, layer_is_incompressible_tuple : sequence of bool
        Per-layer static and incompressible flags; passed through to the output.
    shear_rheology_model_tuple, bulk_rheology_model_tuple : rheology model, str, or sequence of them
        ``TidalPy.rheology_x`` model instances (or model names accepted by ``make_rheology``). A single
        model or name applies to every layer; a sequence supplies one per layer.
    radius_fraction_tuple, thickness_fraction_tuple, volume_fraction_tuple : sequence of float, optional
        Layer sizes as cumulative upper-radius fractions (last entry 1), thickness fractions (sum 1),
        or volume fractions (sum 1). Exactly one must be provided.
    slices_tuple : sequence of int, optional
        Number of radial slices per layer (at least 5 each). Overrides `slice_per_layer`.
    slice_per_layer : int, default 10
        Number of radial slices used for every layer when `slices_tuple` is not provided.
    perform_checks : bool, default True
        Kept for signature compatibility with the classic builder; the native builder always
        validates its inputs.

    Returns
    -------
    PlanetBuildData
        Named tuple whose fields are the positional arguments of `radial_solver`.
    """
    cdef size_t num_layers = <size_t>len(density_tuple)
    if num_layers == 0:
        raise ValueError("At least one layer is required.")

    # Per-layer scalar inputs.
    cdef vector[double] density_vec, bulk_vec, shear_vec, bulk_visc_vec, shear_visc_vec
    _fill_double_vector(_as_layer_tuple(density_tuple, num_layers, "density_tuple"),
                        "density_tuple", density_vec)
    _fill_double_vector(_as_layer_tuple(static_bulk_modulus_tuple, num_layers, "static_bulk_modulus_tuple"),
                        "static_bulk_modulus_tuple", bulk_vec)
    _fill_double_vector(_as_layer_tuple(static_shear_modulus_tuple, num_layers, "static_shear_modulus_tuple"),
                        "static_shear_modulus_tuple", shear_vec)
    _fill_double_vector(_as_layer_tuple(bulk_viscosity_tuple, num_layers, "bulk_viscosity_tuple"),
                        "bulk_viscosity_tuple", bulk_visc_vec)
    _fill_double_vector(_as_layer_tuple(shear_viscosity_tuple, num_layers, "shear_viscosity_tuple"),
                        "shear_viscosity_tuple", shear_visc_vec)

    # Pass-through layer descriptors.
    cdef tuple layer_types = _as_layer_tuple(layer_type_tuple, num_layers, "layer_type_tuple")
    cdef tuple is_static = _as_layer_tuple(layer_is_static_tuple, num_layers, "layer_is_static_tuple")
    cdef tuple is_incompressible = _as_layer_tuple(
        layer_is_incompressible_tuple, num_layers, "layer_is_incompressible_tuple")

    # Layer sizes: exactly one description.
    cdef int num_fraction_inputs = (
        (radius_fraction_tuple is not None) + (thickness_fraction_tuple is not None) +
        (volume_fraction_tuple is not None))
    if num_fraction_inputs != 1:
        raise ValueError(
            "Provide exactly one of `thickness_fraction_tuple`, `radius_fraction_tuple`, or "
            "`volume_fraction_tuple`.")
    cdef vector[double] fraction_vec, thickness_vec
    if thickness_fraction_tuple is not None:
        _fill_double_vector(_as_layer_tuple(thickness_fraction_tuple, num_layers, "thickness_fraction_tuple"),
                            "thickness_fraction_tuple", thickness_vec)
    elif radius_fraction_tuple is not None:
        _fill_double_vector(_as_layer_tuple(radius_fraction_tuple, num_layers, "radius_fraction_tuple"),
                            "radius_fraction_tuple", fraction_vec)
        c_thickness_from_radius_fractions(fraction_vec, thickness_vec)
    else:
        _fill_double_vector(_as_layer_tuple(volume_fraction_tuple, num_layers, "volume_fraction_tuple"),
                            "volume_fraction_tuple", fraction_vec)
        c_thickness_from_volume_fractions(planet_radius, fraction_vec, thickness_vec)

    # Slice counts.
    cdef vector[size_t] slices_vec
    if slices_tuple is not None:
        _fill_size_vector(_as_layer_tuple(slices_tuple, num_layers, "slices_tuple"), "slices_tuple", slices_vec)
    else:
        slices_vec.resize(num_layers, slice_per_layer)

    # Rheology models (one per layer; a single model is broadcast).
    cdef vector[const c_RheologyBase*] shear_rheo_ptrs, bulk_rheo_ptrs
    cdef list shear_keep_alive = cy_resolve_rheology_bylayer(
        shear_rheology_model_tuple, num_layers, "shear_rheology_model_tuple", shear_rheo_ptrs)
    cdef list bulk_keep_alive = cy_resolve_rheology_bylayer(
        bulk_rheology_model_tuple, num_layers, "bulk_rheology_model_tuple", bulk_rheo_ptrs)

    cdef c_RadialSolverInputs inputs
    c_build_rs_input_homogeneous_layers(
        planet_radius, forcing_frequency,
        density_vec, bulk_vec, shear_vec, bulk_visc_vec, shear_visc_vec,
        thickness_vec, slices_vec, shear_rheo_ptrs, bulk_rheo_ptrs, inputs)

    return _build_outputs(inputs, layer_types, is_static, is_incompressible)


def build_rs_input_from_data(
        double forcing_frequency,
        radius_array,
        density_array,
        static_bulk_modulus_array,
        static_shear_modulus_array,
        bulk_viscosity_array,
        shear_viscosity_array,
        layer_upper_radius_tuple,
        layer_type_tuple,
        layer_is_static_tuple,
        layer_is_incompressible_tuple,
        shear_rheology_model_tuple,
        bulk_rheology_model_tuple,
        cpp_bool perform_checks=True,
        cpp_bool warnings=True):
    """Build radial-solver inputs from radially resolved data (for example an external interior model).

    The radius grid is copied and repaired where the solver's requirements are not met: a slice is
    inserted at r = 0 if missing, at each layer base when the previous layer's upper radius is not
    repeated, and at each layer top when it is absent. Inserted slices copy the properties of the
    neighbouring provided slice (a base copies the layer's first slice, a top copies the slice below).
    An interface radius listed only once counts as the top of the lower layer. Each repair is logged as
    a warning when `warnings` is set. The complex moduli are evaluated at `forcing_frequency` with the
    supplied ``rheology_x`` models.

    Parameters
    ----------
    forcing_frequency : float
        Tidal forcing frequency [rad s-1].
    radius_array : array-like of float
        Ascending radius grid [m]; the last entry is the planet radius.
    density_array, static_bulk_modulus_array, static_shear_modulus_array : array-like of float
        Density [kg m-3] and unrelaxed bulk and shear moduli [Pa] at each radius.
    bulk_viscosity_array, shear_viscosity_array : array-like of float
        Bulk and shear viscosities [Pa s] at each radius.
    layer_upper_radius_tuple : sequence of float
        Upper radius of each layer [m], increasing; the last entry must equal the planet radius.
    layer_type_tuple : sequence of str
        Per-layer type ("solid" or "liquid"); passed through to the output.
    layer_is_static_tuple, layer_is_incompressible_tuple : sequence of bool
        Per-layer static and incompressible flags; passed through to the output.
    shear_rheology_model_tuple, bulk_rheology_model_tuple : rheology model, str, or sequence of them
        ``TidalPy.rheology_x`` model instances (or model names accepted by ``make_rheology``). A single
        model or name applies to every layer; a sequence supplies one per layer.
    perform_checks : bool, default True
        Kept for signature compatibility with the classic builder; the native builder always
        validates its inputs.
    warnings : bool, default True
        Log a warning for each grid repair.

    Returns
    -------
    PlanetBuildData
        Named tuple whose fields are the positional arguments of `radial_solver`.
    """
    cdef size_t num_layers = <size_t>len(layer_upper_radius_tuple)
    if num_layers == 0:
        raise ValueError("At least one layer is required.")

    cdef vector[double] radius_vec, density_vec, bulk_vec, shear_vec, bulk_visc_vec, shear_visc_vec
    _fill_vector_from_array(radius_array, "radius_array", radius_vec)
    _fill_vector_from_array(density_array, "density_array", density_vec)
    _fill_vector_from_array(static_bulk_modulus_array, "static_bulk_modulus_array", bulk_vec)
    _fill_vector_from_array(static_shear_modulus_array, "static_shear_modulus_array", shear_vec)
    _fill_vector_from_array(bulk_viscosity_array, "bulk_viscosity_array", bulk_visc_vec)
    _fill_vector_from_array(shear_viscosity_array, "shear_viscosity_array", shear_visc_vec)

    cdef vector[double] upper_radius_vec
    _fill_double_vector(_as_layer_tuple(layer_upper_radius_tuple, num_layers, "layer_upper_radius_tuple"),
                        "layer_upper_radius_tuple", upper_radius_vec)

    cdef tuple layer_types = _as_layer_tuple(layer_type_tuple, num_layers, "layer_type_tuple")
    cdef tuple is_static = _as_layer_tuple(layer_is_static_tuple, num_layers, "layer_is_static_tuple")
    cdef tuple is_incompressible = _as_layer_tuple(
        layer_is_incompressible_tuple, num_layers, "layer_is_incompressible_tuple")

    cdef vector[const c_RheologyBase*] shear_rheo_ptrs, bulk_rheo_ptrs
    cdef list shear_keep_alive = cy_resolve_rheology_bylayer(
        shear_rheology_model_tuple, num_layers, "shear_rheology_model_tuple", shear_rheo_ptrs)
    cdef list bulk_keep_alive = cy_resolve_rheology_bylayer(
        bulk_rheology_model_tuple, num_layers, "bulk_rheology_model_tuple", bulk_rheo_ptrs)

    cdef c_RadialSolverInputs inputs
    c_build_rs_input_from_data(
        forcing_frequency,
        radius_vec, density_vec, bulk_vec, shear_vec, bulk_visc_vec, shear_visc_vec,
        upper_radius_vec, shear_rheo_ptrs, bulk_rheo_ptrs, warnings, inputs)

    return _build_outputs(inputs, layer_types, is_static, is_incompressible)
