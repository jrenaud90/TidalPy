# distutils: language = c++
# cython: boundscheck=False, wraparound=False, nonecheck=False, cdivision=True, initializedcheck=False
"""
prem.pyx
Fast loading and layer auto-detection for PREM-like radial data files.

A PREM-like data file is a delimited table (comma, tab, or whitespace separated)
with columns:

    radius [km], density [kg/m^3], V_p [m/s], V_s [m/s]

and, optionally, two more columns:

    shear viscosity [Pa s], bulk viscosity [Pa s]

From the density and seismic velocities the static shear and bulk moduli are
derived (``mu = rho * Vs^2``; ``K = rho * (Vp^2 - 4/3 Vs^2)``). The radial profile
is then scanned, from the center outward, to split it into layers: a slice with
zero shear modulus is liquid, non-zero is solid, and every solid<->liquid
transition starts a new layer. This scan is the performance-sensitive part of a
PREM planet build and is written in Cython.

All radii are converted to MKS (meters) on load; the file's own ordering
(surface-first or center-first) does not matter because the arrays are sorted
ascending in radius.
"""

from libcpp cimport bool as cpp_bool

import numpy as np
cimport numpy as cnp


# Shear moduli below this threshold [Pa] are treated as zero (liquid). A tiny
# positive floor avoids classifying round-off noise in Vs as a solid.
DEF _DEFAULT_SHEAR_FLOOR_PA = 1.0e-6


def _detect_delimiter(str file_path):
    """Return the delimiter of a data file: ``','``, ``'\\t'``, or None (whitespace)."""
    with open(file_path, "r") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            if "," in line:
                return ","
            if "\t" in line:
                return "\t"
            return None  # whitespace-delimited (np.genfromtxt default)
    return None


def load_prem_arrays(str file_path):
    """Read a PREM-like data file and return its MKS arrays plus derived moduli.

    Parameters
    ----------
    file_path : str
        Path to the delimited data file. The delimiter (comma, tab, or
        whitespace) is auto-detected. Lines beginning with ``#`` are ignored.

    Returns
    -------
    dict
        Keys (all ``numpy.ndarray`` sorted ascending in radius, MKS):
        ``radius_m``, ``density_kg_m3``, ``vp_m_s``, ``vs_m_s``,
        ``shear_modulus_pa`` (``rho*Vs^2``), ``bulk_modulus_pa``
        (``rho*(Vp^2 - 4/3 Vs^2)``), and ``shear_viscosity_pas`` /
        ``bulk_viscosity_pas`` (each ``None`` if the file has no such column).

    Raises
    ------
    ValueError
        If the file has fewer than four columns or an unsupported column count.
    """
    delimiter = _detect_delimiter(file_path)
    data = np.genfromtxt(file_path, delimiter=delimiter, comments="#", dtype=np.float64)
    data = np.atleast_2d(data)
    if data.shape[1] < 4:
        raise ValueError(
            f"PREM-like data file '{file_path}' must have at least 4 columns "
            "(radius_km, density, Vp, Vs); found "
            f"{data.shape[1]}.")
    if data.shape[1] not in (4, 6):
        raise ValueError(
            f"PREM-like data file '{file_path}' must have 4 columns "
            "(radius_km, density, Vp, Vs) or 6 columns (+ shear viscosity, bulk "
            f"viscosity); found {data.shape[1]}.")

    # Sort center-to-surface (ascending radius) so the scan runs bottom-to-top.
    order = np.argsort(data[:, 0], kind="stable")
    data = data[order]

    radius_m      = np.ascontiguousarray(data[:, 0] * 1.0e3)   # km -> m
    density_kg_m3 = np.ascontiguousarray(data[:, 1])
    vp_m_s        = np.ascontiguousarray(data[:, 2])
    vs_m_s        = np.ascontiguousarray(data[:, 3])

    # Static moduli from density + seismic velocities.
    shear_modulus_pa = np.ascontiguousarray(density_kg_m3 * vs_m_s * vs_m_s)
    bulk_modulus_pa  = np.ascontiguousarray(density_kg_m3 * (vp_m_s * vp_m_s - (4.0 / 3.0) * vs_m_s * vs_m_s))

    shear_viscosity_pas = None
    bulk_viscosity_pas  = None
    if data.shape[1] == 6:
        shear_viscosity_pas = np.ascontiguousarray(data[:, 4])
        bulk_viscosity_pas  = np.ascontiguousarray(data[:, 5])

    return {
        "radius_m":            radius_m,
        "density_kg_m3":       density_kg_m3,
        "vp_m_s":              vp_m_s,
        "vs_m_s":              vs_m_s,
        "shear_modulus_pa":    shear_modulus_pa,
        "bulk_modulus_pa":     bulk_modulus_pa,
        "shear_viscosity_pas": shear_viscosity_pas,
        "bulk_viscosity_pas":  bulk_viscosity_pas,
    }


cpdef list detect_layer_boundaries(
        double[::1] radius_m,
        double[::1] shear_modulus_pa,
        double shear_floor_pa=_DEFAULT_SHEAR_FLOOR_PA):
    """Split a radial profile into layers by shear modulus (solid vs liquid).

    Scans the (ascending-radius) profile from the center outward. A slice with
    shear modulus at or below ``shear_floor_pa`` is liquid; above it is solid.
    Each solid<->liquid transition begins a new layer.

    Parameters
    ----------
    radius_m : memoryview of double
        Radii [m], ascending (center to surface).
    shear_modulus_pa : memoryview of double
        Shear modulus [Pa] at each radius (same length as ``radius_m``).
    shear_floor_pa : float, optional
        Shear moduli at or below this are treated as zero (liquid).

    Returns
    -------
    list of tuple
        One ``(start_index, end_index, is_solid)`` per layer, inner to outer.
        ``start_index`` / ``end_index`` are inclusive indices into the arrays;
        ``is_solid`` is a Python bool.

    Notes
    -----
    Profiles like PREM mark a phase boundary with two points at the *same* radius
    (e.g. the inner-core/outer-core boundary), which would otherwise produce
    zero-thickness "layers". Such degenerate runs (whose outer radius does not exceed
    their inner radius) are absorbed into the previous real layer, so only layers
    spanning a non-zero radius interval are returned.
    """
    cdef Py_ssize_t n = radius_m.shape[0]
    cdef list raw = []
    cdef list merged = []
    cdef Py_ssize_t i
    cdef Py_ssize_t start = 0
    cdef Py_ssize_t idx_a, idx_b
    cdef cpp_bool current_solid
    cdef cpp_bool slice_solid

    if n == 0:
        return []

    # First pass: maximal runs of constant solidity, recorded as Python tuples
    # (start_index, end_index, is_solid).
    current_solid = shear_modulus_pa[0] > shear_floor_pa
    for i in range(1, n):
        slice_solid = shear_modulus_pa[i] > shear_floor_pa
        if slice_solid != current_solid:
            raw.append((start, i - 1, bool(current_solid)))
            start = i
            current_solid = slice_solid
    raw.append((start, n - 1, bool(current_solid)))

    # Second pass: absorb zero-thickness (duplicate-radius) runs into the previous
    # real layer, keeping that layer's solidity. A reference to the current layer is
    # tracked explicitly (a negative list index is undefined under wraparound=False).
    last_layer = None
    for run in raw:
        idx_a = run[0]
        idx_b = run[1]
        if last_layer is not None and radius_m[idx_b] <= radius_m[idx_a]:
            last_layer[1] = run[1]
        else:
            last_layer = [run[0], run[1], run[2]]
            merged.append(last_layer)

    return [(item[0], item[1], bool(item[2])) for item in merged]
