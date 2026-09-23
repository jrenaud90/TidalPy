"""Radial data files: a PREM-like profile in, a purely interpolated planet out.

A world config may describe its interior with a radial profile rather than by listing each layer's
geometry and material. The profile is either a delimited file

    data_file = "PREM.csv"

or, from Python, a mapping of arrays

    build_world({..., "data": {"radius_km": [...], "density": [...], "vp": [...], "vs": [...]}})

A profile says how many layers the world has, where their boundaries are, and what each is made of. It
does not say everything a layer can carry: it holds no complex moduli, so a rheology is still named in a
``[layers.<name>]`` table, as are a cooling model and radiogenics. Such a table refines one detected
layer, and only the layers being refined need one.

Both go through :func:`load_radial_data`, which normalizes them to MKS arrays ascending in radius. A
profile must carry a radius (or a depth), a density, and the seismic velocities V_p and V_s, from which
the static moduli follow

    mu = rho V_s^2                K = rho (V_p^2 - 4/3 V_s^2)

A profile may instead give those moduli directly, and may add a shear and a bulk viscosity. Without
viscosities the layers are elastic: nothing dissipates, and no viscosity or partial-melt model is built.

Columns are found by name, so their order does not matter and they may state their units; a name whose
unit is not one this reader converts is taken to be MKS already. A file with no header is read
positionally in the canonical order. :func:`detect_layer_boundaries` then splits the profile at every
solid/liquid transition, and the world builder hands each layer's slice to an interpolated EOS, which
owns those arrays. That EOS is the only place a radial grid persists.
"""

import os
import re
from typing import Optional, Union

import numpy as np

# Shear moduli at or below this threshold [Pa] are treated as zero (liquid). A tiny positive floor
# avoids classifying round-off noise in V_s as a solid.
DEFAULT_SHEAR_FLOOR_PA = 1.0e-6

# The canonical MKS quantities of a radial profile. Each maps to the spellings accepted for its
# column, reduced to lowercase words joined by underscores (so "V_p", "Vp [m/s]", and "VP" all reduce
# to "vp"). A name may append its unit to any of these.
_COLUMN_ALIASES = {
    "radius":          ("radius", "r", "rad", "radii"),
    "density":         ("density", "rho", "dens"),
    "vp":              ("vp", "v_p", "vpv", "p_velocity", "p_wave_velocity", "compressional_velocity"),
    "vs":              ("vs", "v_s", "vsv", "s_velocity", "s_wave_velocity", "shear_velocity"),
    "shear_viscosity": ("shear_viscosity", "viscosity", "eta", "eta_shear", "shear_eta", "visc"),
    "bulk_viscosity":  ("bulk_viscosity", "eta_bulk", "bulk_eta", "zeta"),
    # Alternatives to the velocities: the static moduli themselves.
    "shear_modulus":   ("shear_modulus", "mu", "shear", "rigidity"),
    "bulk_modulus":    ("bulk_modulus", "k", "bulk", "incompressibility"),
    # An alternative to the radius, converted with the world's surface radius.
    "depth":           ("depth", "z"),
}

# The order a headerless file must list its columns in. The last two are optional.
POSITIONAL_ORDER = ("radius", "density", "vp", "vs", "shear_viscosity", "bulk_viscosity")
_REQUIRED_POSITIONAL = 4

# Unit suffixes this reader converts, by the kind of quantity carrying them. Any other suffix (``pa``,
# ``kg_m3``, ``pas``) names a unit that is already MKS.
_LENGTH_UNITS = {"m": 1.0, "meter": 1.0, "meters": 1.0, "metre": 1.0, "metres": 1.0,
                 "km": 1.0e3, "kilometer": 1.0e3, "kilometers": 1.0e3,
                 "kilometre": 1.0e3, "kilometres": 1.0e3}
_VELOCITY_UNITS = {"m_s": 1.0, "ms": 1.0, "m_s1": 1.0, "m_per_s": 1.0,
                   "km_s": 1.0e3, "kms": 1.0e3, "km_s1": 1.0e3, "km_per_s": 1.0e3}

# A radius or depth given without a unit is read as kilometers below this value [m] and as meters at
# or above it. A planet large enough for this library is at least 100 km in radius, and no radius in
# kilometers reaches 100000, so the two ranges cannot overlap.
_KM_CUTOFF_M = 1.0e5


def _normalize_name(name: str) -> str:
    """Reduce a column name to lowercase words joined by single underscores (``"Vp [km/s]"`` -> ``"vp_km_s"``)."""
    return re.sub(r"[^0-9a-zA-Z]+", "_", str(name).strip().lower()).strip("_")


def _match_quantity(name: str):
    """Map a column name to ``(quantity, unit)``, or ``(None, None)`` when it names nothing this reader reads.

    An exact alias wins; otherwise the longest alias the name starts with claims it, and the rest of
    the name is its unit. A unit this reader does not convert (``pa``, ``kg_m3``) yields ``None``,
    meaning the column is MKS already.
    """
    stem = _normalize_name(name)
    best = None
    for quantity, aliases in _COLUMN_ALIASES.items():
        for alias in aliases:
            if stem == alias:
                return quantity, None
            if stem.startswith(alias + "_") and (best is None or len(alias) > len(best[1])):
                best = (quantity, alias)
    if best is None:
        return None, None
    quantity, alias = best
    unit = stem[len(alias) + 1:]
    if quantity in ("radius", "depth"):
        return quantity, _LENGTH_UNITS.get(unit)
    if quantity in ("vp", "vs"):
        return quantity, _VELOCITY_UNITS.get(unit)
    return quantity, None


def _length_to_meters(values: np.ndarray, factor: Optional[float]) -> np.ndarray:
    """Convert a radius or depth column to meters, from its stated unit or, failing that, its magnitude."""
    if factor is not None:
        return values * factor
    largest = float(np.max(np.abs(values))) if values.size else 0.0
    if largest < _KM_CUTOFF_M:
        return values * 1.0e3   # kilometers: the usual spelling of a PREM-like table
    return values


# =====================================================================================================================
# Reading a delimited file
# =====================================================================================================================
def _detect_delimiter(line: str) -> Optional[str]:
    """Return the delimiter of a data line: ``','``, ``';'``, ``'\\t'``, or None for whitespace."""
    for candidate in (",", ";", "\t"):
        if candidate in line:
            return candidate
    return None


def _read_table(file_path: str):
    """Read a delimited data file, returning ``(data, header)``: a 2-D float array and its column names or None.

    The delimiter (comma, semicolon, tab, or whitespace) is taken from the first non-comment line.
    A header is a leading non-numeric row or, failing that, the last comment line before the data
    that has one field per column and names a quantity in every one of them.
    """
    comments = []
    header = None
    delimiter = None
    rows = []
    with open(file_path, "r") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("#"):
                if not rows:
                    comments.append(line.lstrip("#").strip())
                continue
            if delimiter is None and not rows and header is None:
                delimiter = _detect_delimiter(line)
            fields = [field.strip() for field in (line.split(delimiter) if delimiter else line.split())]
            fields = [field for field in fields if field]
            if not fields:
                continue
            try:
                rows.append([float(field) for field in fields])
            except ValueError:
                if header is None and not rows:
                    header = fields         # a leading non-numeric row names the columns
                    continue
                raise ValueError(
                    f"Radial data file '{file_path}' line {line_number} is not numeric: {line!r}.") from None

    if not rows:
        raise ValueError(f"Radial data file '{file_path}' holds no data rows.")
    widths = {len(row) for row in rows}
    if len(widths) > 1:
        raise ValueError(
            f"Radial data file '{file_path}' has rows of differing widths ({sorted(widths)}); every "
            "row must have the same number of columns.")
    data = np.array(rows, dtype=np.float64)

    if header is None:
        # A comment line is only a guess at a header, so take it as one only when it has a field per
        # column and every one of them names a quantity. Prose about the data often has the right
        # number of commas and a word or two in common with it, which is not the same thing.
        for comment in reversed(comments):
            fields = [field.strip() for field in (comment.split(delimiter) if delimiter else comment.split())]
            fields = [field for field in fields if field]
            if len(fields) == data.shape[1] and all(_match_quantity(field)[0] for field in fields):
                header = fields
                break
    if header is not None and len(header) != data.shape[1]:
        header = None
    return data, header


def _columns_from_file(file_path: str) -> dict:
    """Read a delimited file into ``{quantity: (values, unit_factor)}``, by column name or by position."""
    data, header = _read_table(file_path)
    columns = {}
    if header is not None:
        unmatched = []
        for index, name in enumerate(header):
            quantity, factor = _match_quantity(name)
            if quantity is None:
                unmatched.append(name)
                continue
            if quantity in columns:
                raise ValueError(
                    f"Radial data file '{file_path}' names the {quantity} column more than once "
                    f"(one of them is '{name}').")
            columns[quantity] = (data[:, index], factor)
        if not columns:
            raise ValueError(
                f"Radial data file '{file_path}' has a header naming nothing this reader knows "
                f"({', '.join(unmatched)}). Expected names such as radius_km, density, vp, vs.")
        return columns

    if data.shape[1] < _REQUIRED_POSITIONAL:
        raise ValueError(
            f"Radial data file '{file_path}' has no header, so its columns are read in the order "
            f"{', '.join(POSITIONAL_ORDER[:_REQUIRED_POSITIONAL])} (then the optional shear and bulk "
            f"viscosity); that needs at least {_REQUIRED_POSITIONAL} columns but the file has "
            f"{data.shape[1]}.")
    if data.shape[1] > len(POSITIONAL_ORDER):
        raise ValueError(
            f"Radial data file '{file_path}' has no header and {data.shape[1]} columns; a headerless "
            f"file may have at most {len(POSITIONAL_ORDER)} ({', '.join(POSITIONAL_ORDER)}). Add a "
            "header line naming the columns.")
    return {POSITIONAL_ORDER[index]: (data[:, index], None) for index in range(data.shape[1])}


def _columns_from_mapping(mapping) -> dict:
    """Read a mapping of arrays into ``{quantity: (values, unit_factor)}``, keyed as a file header is."""
    columns = {}
    unmatched = []
    for name, values in mapping.items():
        quantity, factor = _match_quantity(name)
        if quantity is None:
            unmatched.append(str(name))
            continue
        if quantity in columns:
            raise ValueError(f"Radial data names the {quantity} column more than once (one of them is '{name}').")
        array = np.ascontiguousarray(values, dtype=np.float64).ravel()
        columns[quantity] = (array, factor)
    if not columns:
        raise ValueError(
            "Radial data names nothing this reader knows "
            f"({', '.join(unmatched) if unmatched else 'the mapping is empty'}). Expected keys such as "
            "radius_km, density, vp, vs.")
    return columns


# =====================================================================================================================
# Normalization
# =====================================================================================================================
def load_radial_data(source: Union[str, dict], surface_radius: Optional[float] = None) -> dict:
    """Read a radial profile from a file or a mapping of arrays and return its MKS arrays.

    Parameters
    ----------
    source : str or dict
        Path to a delimited data file, or a mapping of column name to sequence. Names are matched
        case- and punctuation-insensitively against the known spellings and may carry a unit
        (``radius_km``, ``vp_km_s``). A file with no header is read positionally as radius, density,
        V_p, V_s, shear viscosity, bulk viscosity.
    surface_radius : float, optional
        The world's radius [m]. Needed only to convert a profile given as a depth.

    Returns
    -------
    dict
        ``radius_m``, ``density_kg_m3``, ``shear_modulus_pa``, and ``bulk_modulus_pa`` always;
        ``vp_m_s`` / ``vs_m_s`` when the profile gave velocities and ``shear_viscosity_pas`` /
        ``bulk_viscosity_pas`` when it gave viscosities (``None`` otherwise). Every array is
        contiguous ``float64``, sorted ascending in radius.

    Raises
    ------
    ValueError
        A required quantity is missing, the columns disagree in length, or a value is unphysical.
    """
    if isinstance(source, (str, bytes, os.PathLike)):
        columns = _columns_from_file(os.fspath(source))
        where = f" in '{source}'"
    elif hasattr(source, "items"):
        columns = _columns_from_mapping(source)
        where = ""
    else:
        raise TypeError(
            "Radial data must be a path to a data file or a mapping of column name to array; got "
            f"{type(source).__name__}.")

    lengths = {quantity: values.shape[0] for quantity, (values, _) in columns.items()}
    if len(set(lengths.values())) > 1:
        raise ValueError(f"Radial data{where} has columns of differing lengths: {lengths}.")

    # ---- the radius ------------------------------------------------------------------------------------------
    if "radius" in columns:
        radius = _length_to_meters(*columns["radius"])
    elif "depth" in columns:
        if surface_radius is None:
            raise ValueError(
                f"Radial data{where} gives a depth rather than a radius, so the world's 'radius_m' is "
                "needed to convert it; set that key on the world.")
        radius = float(surface_radius) - _length_to_meters(*columns["depth"])
    else:
        raise ValueError(f"Radial data{where} has no radius (or depth) column.")

    if "density" not in columns:
        raise ValueError(f"Radial data{where} has no density column.")
    density = columns["density"][0]

    # ---- the static moduli, from the velocities or given outright --------------------------------------------
    vp = None
    vs = None
    if "vp" in columns and "vs" in columns:
        vp = _velocity(columns["vp"])
        vs = _velocity(columns["vs"])
        shear_modulus = density * vs * vs
        bulk_modulus  = density * (vp * vp - (4.0 / 3.0) * vs * vs)
    elif "shear_modulus" in columns and "bulk_modulus" in columns:
        shear_modulus = columns["shear_modulus"][0]
        bulk_modulus  = columns["bulk_modulus"][0]
    else:
        raise ValueError(
            f"Radial data{where} must give the seismic velocities (columns vp and vs) or the static "
            "moduli themselves (columns shear_modulus and bulk_modulus); it gives neither pair.")

    arrays = {
        "radius_m":            radius,
        "density_kg_m3":       density,
        "vp_m_s":              vp,
        "vs_m_s":              vs,
        "shear_modulus_pa":    shear_modulus,
        "bulk_modulus_pa":     bulk_modulus,
        "shear_viscosity_pas": columns["shear_viscosity"][0] if "shear_viscosity" in columns else None,
        "bulk_viscosity_pas":  columns["bulk_viscosity"][0] if "bulk_viscosity" in columns else None,
    }
    if radius.size < 2:
        raise ValueError(f"Radial data{where} has {radius.size} row(s); a profile needs at least 2.")

    # Sort center-to-surface so the layer scan runs bottom-to-top. A surface-first profile is reversed
    # first: the stable sort then keeps the lower layer's row ahead of the upper layer's at a duplicated
    # boundary radius, so each layer's arrays end (and start) with its own values.
    if arrays["radius_m"][0] > arrays["radius_m"][-1]:
        arrays = {key: (None if value is None else value[::-1]) for key, value in arrays.items()}
    order = np.argsort(arrays["radius_m"], kind="stable")
    arrays = {key: (None if value is None else np.ascontiguousarray(value[order]))
              for key, value in arrays.items()}

    _validate_profile(arrays, where)
    return arrays


def _velocity(column) -> np.ndarray:
    """Convert a seismic velocity column to m/s; one with no stated unit is already MKS."""
    values, factor = column
    return values if factor is None else values * factor


def _validate_profile(arrays: dict, where: str) -> None:
    """Reject a profile that cannot describe a planet (non-finite, negative, or otherwise unphysical)."""
    for key, values in arrays.items():
        if values is not None and not np.all(np.isfinite(values)):
            raise ValueError(f"Radial data{where} has non-finite values in its {key} column.")
    if arrays["radius_m"][0] < 0.0:
        raise ValueError(f"Radial data{where} has a negative radius.")
    if arrays["radius_m"][-1] <= arrays["radius_m"][0]:
        raise ValueError(
            f"Radial data{where} spans no radius: every row sits at "
            f"{float(arrays['radius_m'][0]):.6g} m.")
    if np.any(arrays["density_kg_m3"] <= 0.0):
        raise ValueError(f"Radial data{where} has a non-positive density.")
    if np.any(arrays["shear_modulus_pa"] < 0.0):
        raise ValueError(f"Radial data{where} has a negative shear modulus.")
    if np.any(arrays["bulk_modulus_pa"] <= 0.0):
        raise ValueError(f"Radial data{where} has a non-positive bulk modulus.")
    for key in ("shear_viscosity_pas", "bulk_viscosity_pas"):
        if arrays[key] is not None and np.any(arrays[key] <= 0.0):
            raise ValueError(f"Radial data{where} has a non-positive {key[:-4].replace('_', ' ')}.")


# =====================================================================================================================
# Layer detection
# =====================================================================================================================
def detect_layer_boundaries(radius, shear_modulus, shear_floor: float = DEFAULT_SHEAR_FLOOR_PA) -> list:
    """Split a radial profile into layers by shear modulus (solid vs liquid).

    Scans the (ascending-radius) profile from the center outward. A slice with shear modulus at or
    below ``shear_floor`` is liquid; above it is solid. Each solid<->liquid transition begins a new
    layer.

    Parameters
    ----------
    radius : array-like of float
        Radii [m], ascending (center to surface).
    shear_modulus : array-like of float
        Shear modulus [Pa] at each radius (same length as ``radius``).
    shear_floor : float, optional
        Shear moduli at or below this are treated as zero (liquid).

    Returns
    -------
    list of tuple
        One ``(start_index, end_index, is_solid)`` per layer, inner to outer. The indices are
        inclusive; ``is_solid`` is a Python bool.

    Notes
    -----
    Profiles like PREM mark a phase boundary with two points at the *same* radius (the
    inner-core/outer-core boundary, for one), which would otherwise produce zero-thickness "layers".
    Such degenerate runs (whose outer radius does not exceed their inner radius) are absorbed into
    the previous real layer, so every returned layer spans a non-zero radius interval.
    """
    radius = np.ascontiguousarray(radius, dtype=np.float64)
    shear_modulus = np.ascontiguousarray(shear_modulus, dtype=np.float64)
    if radius.size == 0:
        return []

    # Maximal runs of constant solidity: a run starts wherever solidity changes.
    solid = shear_modulus > shear_floor
    starts = np.concatenate(([0], np.flatnonzero(np.diff(solid)) + 1))
    ends = np.concatenate((starts[1:] - 1, [radius.size - 1]))

    # Absorb zero-thickness (duplicate-radius) runs into the previous real layer, keeping that layer's
    # solidity.
    layers = []
    for start, end in zip(starts, ends):
        if layers and radius[end] <= radius[start]:
            layers[-1][1] = int(end)
        else:
            layers.append([int(start), int(end), bool(solid[start])])

    # A run at the very center has no previous layer to be absorbed into, so it is absorbed forward
    # into the first real one instead (a profile whose first two rows sit at radius 0).
    while len(layers) > 1 and radius[layers[0][1]] <= radius[layers[0][0]]:
        layers[1][0] = layers[0][0]
        layers.pop(0)
    return [(start, end, is_solid) for start, end, is_solid in layers]
