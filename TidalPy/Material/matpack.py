"""MatPack: the named materials TidalPy ships, and the preset tables that start from them.

Each file in ``TidalPy/MatPack`` holds one material config table (the form ``Material.get_config_dict()`` returns and
``make_material`` takes) and three metadata keys: ``schema_version``, ``description``, and ``category`` (one of
:data:`CATEGORIES`). Its references are comments in the file. The files are installed copy-if-absent into
``<documents>/TidalPy/<version>/Materials`` (see :class:`TidalPy.Utilities.data_pack.DataPack`), where they can be
edited; an edited copy is preferred to the packaged file.

A material is given in one of three forms:

- a name, ``"peridotite"``;
- a preset with overrides, ``{"preset": "peridotite", "latent_heat_j_kg": 3.0e5, "solid": {...}}``;
- a full material table without ``preset``, used as given.

Overrides merge into the preset's table one table at a time, so they name only what they change. A model table
(``eos``, ``shear_viscosity``, ``solidus``, ...) that names a different model than the preset's replaces the preset's
table instead, since another model reads other keys; ``None`` removes a slot (from Python only, TOML having no null).
A ``solid``, ``liquid``, or ``melting`` table may itself name a preset, ``{"preset": "water"}``, and then starts from
that material's table in the same slot. MatPack files use both forms, so a shared phase such as liquid water, or a
material's melting curves, is written once.
"""

import copy
import difflib
import os
import warnings

import toml

import TidalPy
from TidalPy.configurations import validate_schema_version
from TidalPy.paths import get_materials_dir as _paths_get_materials_dir
from TidalPy.Utilities.classes import canonical_parameter_keys
from TidalPy.Utilities.classes.families import get_family
from TidalPy.Utilities.data_pack import DataPack, user_stacklevel
from TidalPy.Material.material import Material, Phase

# The packaged MatPack directory (read-only source of the materials), relative to the installed package root.
PACKAGED_MATPACK_DIR = os.path.join(os.path.dirname(os.path.abspath(TidalPy.__file__)), "MatPack")

# The key that names a preset, in a material table or in its solid or liquid table.
PRESET_KEY = "preset"

# Keys a MatPack file carries besides its material table.
METADATA_KEYS = ("schema_version", "description", "category")

# The categories a MatPack file names, in the order the documentation lists them.
CATEGORIES = ("simplified", "rocky", "icy", "giant")

# The phase slots of a material table.
_PHASE_SLOTS = ("solid", "liquid")

# The tables of a material table that may name a preset of their own: its phases and its melting table.
_PRESET_SLOTS = _PHASE_SLOTS + ("melting",)

# Model slot -> the family of the model it holds, so an override naming the preset's model by an alias merges into
# the preset's table while one naming another model replaces it.
_SLOT_FAMILIES = {
    "eos":                   "equation of state",
    "shear_modulus":         "shear modulus",
    "shear_viscosity":       "viscosity",
    "bulk_viscosity":        "viscosity",
    "shear_rheology":        "rheology",
    "bulk_rheology":         "rheology",
    "solidus":               "melting curve",
    "liquidus":              "melting curve",
    "weakening":             "melt weakening",
    "bulk_modulus_mixing":   "bulk-modulus mixing",
    "bulk_viscosity_mixing": "bulk-viscosity mixing",
}

# Parsed material files by path, with the text each was parsed from, so a material loaded in a loop does not parse its
# file each time. The text is compared, so an edit is always read.
_PARSED_FILES = {}


def get_materials_dir():
    """Return the user-editable data directory for MatPack materials.

    Thin indirection over :func:`TidalPy.paths.get_materials_dir` so tests can redirect the data directory by patching
    this module attribute.

    Returns
    -------
    str or None
        The ``Materials`` data directory (created if absent); None when it cannot be created, in which case the
        packaged materials are used directly.
    """
    return _paths_get_materials_dir()


# The pack itself. Its getter looks get_materials_dir up on every call, so redirecting that module attribute
# redirects the pack.
MAT_PACK = DataPack(
    "MatPack",
    PACKAGED_MATPACK_DIR,
    lambda: get_materials_dir(),
    (".toml",),
    "stale_matpack_copy",
    "TidalPy.Material.install_matpack(force=True)")


# =====================================================================================================================
# Pack files
# =====================================================================================================================
def install_matpack(force: bool = False):
    """Copy the packaged materials into the data directory (copy-if-absent).

    Parameters
    ----------
    force : bool, optional
        If True, overwrite every data-directory copy with the packaged file (discarding user edits). Default False.

    Returns
    -------
    str or None
        The data directory; None when there is no usable data directory.
    """
    return MAT_PACK.install(force)


def available_materials(category: str = None) -> list:
    """The sorted names of the MatPack materials, optionally of one category.

    Parameters
    ----------
    category : str, optional
        One of :data:`CATEGORIES`; all materials when absent.

    Returns
    -------
    list of str
        Material names, usable wherever a material is named. With a category, a file that cannot be read is left out
        with a warning; it reports its own error when it is loaded.

    Raises
    ------
    ValueError
        For an unknown category.
    """
    if category is not None and category not in CATEGORIES:
        raise ValueError(f"TidalPy: unknown MatPack category '{category}'. Accepted: {', '.join(CATEGORIES)}.")
    names = sorted(MAT_PACK.files(".toml"))
    if category is None:
        return names
    in_category = []
    for name in names:
        try:
            info = material_info(name)
        except (ValueError, OSError, UnicodeDecodeError) as error:
            warnings.warn(f"TidalPy: the MatPack material '{name}' could not be read and is not listed: {error}",
                          stacklevel=user_stacklevel())
            continue
        if info["category"] == category:
            in_category.append(name)
    return in_category


def material_info(name: str) -> dict:
    """A MatPack material's metadata: ``name``, ``description``, ``category``, and the ``path`` it is read from.

    Raises
    ------
    ValueError
        No MatPack material has that name.
    """
    path = p_material_path(name)
    table = p_read_material_file(path)
    return {
        "name": name.lower(),
        "description": table.get("description", ""),
        "category": table.get("category", ""),
        "path": path,
    }


def p_material_path(name: str) -> str:
    """The file a MatPack material is read from: the data-directory copy, else the packaged file."""
    if not isinstance(name, str):
        raise TypeError(f"TidalPy: a material is named by a string, not {type(name).__name__}.")
    if ("/" in name) or ("\\" in name) or (name in ("", ".", "..")):
        raise ValueError(f"TidalPy: '{name}' is not a MatPack material name; a name holds no path.")
    path = MAT_PACK.find(name + ".toml")
    if path is not None:
        return path
    names = sorted(MAT_PACK.files(".toml"))
    close = difflib.get_close_matches(name.lower(), names, n=1)
    hint = f" (did you mean '{close[0]}'?)" if close else ""
    raise ValueError(f"TidalPy: no MatPack material named '{name}'{hint}. Available: {', '.join(names)}.")


def p_read_material_file(path: str) -> dict:
    """A material file's table, metadata included, after its schema-version check; a copy the caller may edit."""
    with open(path, "r", encoding="utf-8") as file:
        text = file.read()
    key = os.path.normcase(os.path.abspath(path))
    cached = _PARSED_FILES.get(key)
    if cached is None or cached[0] != text:
        try:
            parsed = toml.loads(text)
        except toml.TomlDecodeError as error:
            raise ValueError(f"TidalPy: could not parse the material file '{path}': {error}") from error
        validate_schema_version(parsed)
        cached = (text, parsed)
        _PARSED_FILES[key] = cached
    return copy.deepcopy(cached[1])


# =====================================================================================================================
# Presets and overrides
# =====================================================================================================================
def merge_material_tables(base: dict, overrides: dict, kind: str = "material") -> dict:
    """``base`` with ``overrides`` merged over it, leaving both untouched.

    Tables merge key by key, after each parameter in either table is renamed to its config key (a parameter may be
    given by its argument name or its config key), so an override replaces the base's value under either spelling.
    A model table naming a different model than the base's, and a ``solid`` or ``liquid`` table naming a preset,
    replace the base's table, and a value of ``None`` removes the key. Anything else (a number, a list) replaces the
    base value.

    Parameters
    ----------
    base : dict
        A material, phase, or melting table.
    overrides : dict
        The values that win.
    kind : str, optional
        What the tables are: ``"material"`` (default), ``"phase"``, or ``"melting"``.

    Returns
    -------
    dict
        A new, merged table.
    """
    merged = p_canonical_table(copy.deepcopy(base), kind)
    for key, value in p_canonical_table(overrides, kind).items():
        if value is None:
            merged.pop(key, None)
            continue
        base_value = merged.get(key)
        names_preset = (key in _PRESET_SLOTS) and isinstance(value, dict) and (PRESET_KEY in value)
        merges = (isinstance(value, dict) and isinstance(base_value, dict) and not names_preset
                  and not p_names_other_model(key, base_value, value))
        if not merges:
            merged[key] = copy.deepcopy(value)
        elif key in _SLOT_FAMILIES:
            merged[key] = p_merge_model_tables(key, base_value, value)
        else:
            merged[key] = merge_material_tables(base_value, value, p_child_kind(kind, key))
    return merged


def p_child_kind(kind: str, key: str) -> str:
    """What a table inside a ``kind`` table is: a phase, a melting table, or (for anything else) a material."""
    if (kind == "material") and (key in _PHASE_SLOTS):
        return "phase"
    if (kind == "material") and (key == "melting"):
        return "melting"
    return kind


def p_canonical_table(table: dict, kind: str) -> dict:
    """A material, phase, or melting table with its own parameters under their config keys (its tables untouched)."""
    if kind == "material":
        return canonical_parameter_keys(Material, table)
    if kind == "phase":
        return canonical_parameter_keys(Phase, table)
    return dict(table)


def p_model_class(slot: str, model_name: str):
    """The Python class of the model a slot's table names, or None when the name is unknown (building the material
    names that error)."""
    family = get_family(_SLOT_FAMILIES[slot])
    try:
        return family.classes[family.canonical_name(model_name)]
    except (ValueError, KeyError):
        return None


def p_merge_model_tables(slot: str, base_table: dict, override_table: dict) -> dict:
    """Two tables of one model (the override names no other) merged key by key, both keyed by config key."""
    model_name = override_table.get("model", base_table.get("model"))
    model_class = p_model_class(slot, str(model_name)) if model_name is not None else None
    if model_class is not None:
        base_table = canonical_parameter_keys(model_class, base_table)
        override_table = canonical_parameter_keys(model_class, override_table)
    merged = copy.deepcopy(base_table)
    for key, value in override_table.items():
        if value is None:
            merged.pop(key, None)
        else:
            merged[key] = copy.deepcopy(value)
    return merged


def p_names_other_model(slot: str, base_table: dict, override_table: dict) -> bool:
    """Whether an override of a model table names a different model than the table it overrides."""
    if ("model" not in override_table) or ("model" not in base_table):
        return False
    base_model = str(base_table["model"])
    override_model = str(override_table["model"])
    family_name = _SLOT_FAMILIES.get(slot)
    if family_name is None:
        return base_model.lower() != override_model.lower()
    try:
        return not get_family(family_name).same_model(base_model, override_model)
    except ValueError:
        # An unknown model name: the override replaces the table, and building the material names the error.
        return True


def p_resolve(table: dict, chain: tuple) -> dict:
    """A material table with its presets resolved (its own, then each phase's) and its metadata dropped.

    ``chain`` holds the presets being resolved, so a file that names itself, directly or through others, is an
    error rather than an endless loop.
    """
    if not isinstance(table, dict):
        raise TypeError(f"TidalPy: a material is a name or a table, not {type(table).__name__}.")
    table = {key: value for key, value in table.items() if key not in METADATA_KEYS}
    preset = table.pop(PRESET_KEY, None)
    if preset is not None:
        base = p_resolve_preset(preset, chain)
        table = merge_material_tables(base, table)
    for slot in _PRESET_SLOTS:
        slot_table = table.get(slot)
        if isinstance(slot_table, dict) and PRESET_KEY in slot_table:
            table[slot] = p_resolve_slot(slot, slot_table, chain)
    p_check_preset_places(table, "the material table")
    return table


def p_check_preset_places(table: dict, where: str) -> None:
    """Raise for a ``preset`` key below the material and phase levels, where nothing would resolve it."""
    for key, value in table.items():
        if not isinstance(value, dict):
            continue
        if PRESET_KEY in value and key not in _PRESET_SLOTS:
            raise ValueError(f"TidalPy: the '{key}' table of {where} names a '{PRESET_KEY}'; only a material table and "
                             "its 'solid', 'liquid', and 'melting' tables may.")
        p_check_preset_places(value, f"the '{key}' table")


def p_resolve_preset(preset, chain: tuple) -> dict:
    """A MatPack material's resolved table, guarding against a cycle of presets."""
    if not isinstance(preset, str):
        raise TypeError(f"TidalPy: '{PRESET_KEY}' names a MatPack material by a string, not {type(preset).__name__}.")
    name = preset.lower()
    if name in chain:
        raise ValueError(f"TidalPy: the MatPack presets {' -> '.join(chain + (name,))} form a cycle.")
    return p_resolve(p_read_material_file(p_material_path(name)), chain + (name,))


def p_resolve_slot(slot: str, slot_table: dict, chain: tuple) -> dict:
    """A phase or melting table that names a preset: that material's table in the same slot, with the table's other
    keys merged over it."""
    overrides = dict(slot_table)
    preset = overrides.pop(PRESET_KEY)
    material_table = p_resolve_preset(preset, chain)
    if not isinstance(material_table.get(slot), dict):
        noun = "phase" if slot in _PHASE_SLOTS else "table"
        raise ValueError(f"TidalPy: the MatPack material '{preset}' has no '{slot}' {noun} for a '{slot}' table to "
                         "start from.")
    return merge_material_tables(material_table[slot], overrides, "melting" if slot == "melting" else "phase")


def material_config(source, overrides: dict = None) -> dict:
    """The full config table of a material: presets resolved, overrides merged, metadata dropped.

    Parameters
    ----------
    source : str or dict
        A MatPack name, or a material table (with or without a ``preset``).
    overrides : dict, optional
        Merged over the resolved table (see :func:`merge_material_tables`).

    Returns
    -------
    dict
        The table ``make_material`` takes.

    Raises
    ------
    ValueError
        An unknown MatPack name, a cycle of presets, a phase preset with no phase for its slot, a ``preset`` in a
        table that cannot take one, or a result with neither a solid nor a liquid phase.
    TypeError
        A source that is neither a name nor a table.
    """
    if isinstance(source, str):
        table = {PRESET_KEY: source}
    elif isinstance(source, dict):
        table = copy.deepcopy(source)
    else:
        raise TypeError(f"TidalPy: a material is a name or a table, not {type(source).__name__}.")
    resolved = p_resolve(table, ())
    if overrides:
        if PRESET_KEY in overrides:
            raise ValueError(f"TidalPy: overrides cannot name a '{PRESET_KEY}'; name it in the source instead.")
        resolved = p_resolve(merge_material_tables(resolved, overrides), ())
    if not any(isinstance(resolved.get(slot), dict) for slot in _PHASE_SLOTS):
        if "model" in resolved:
            raise ValueError(
                "TidalPy: the material table names a 'model', but a material is made of phases: put the equation of "
                "state in 'solid.eos' (or 'liquid.eos') and the other laws beside it in the phase table.")
        raise ValueError("TidalPy: the material has neither a 'solid' nor a 'liquid' phase.")
    return resolved


def load_material(source, **overrides) -> Material:
    """A material from MatPack, a preset table, or a full table.

    Parameters
    ----------
    source : str or dict
        A MatPack name (``"peridotite"``), or a material table (with or without a ``preset``).
    **overrides
        Merged over the resolved table by config key or slot (``latent_heat_j_kg=3.0e5``,
        ``solid={"shear_viscosity": {...}}``; ``liquid=None, melting=None`` removes the melt).

    Returns
    -------
    Material
        The material, immutable like every model; ``with_parameters`` and ``replace`` give changed copies.

    Raises
    ------
    ValueError
        An unknown name, a preset problem (see :func:`material_config`), or a table the material rejects; the message
        names the source.
    """
    where = f"MatPack material '{source}'" if isinstance(source, str) else "material table"
    try:
        return Material(config=material_config(source, overrides))
    except (ValueError, TypeError) as error:
        detail = str(error)
        if detail.startswith("TidalPy: "):
            detail = detail[len("TidalPy: "):]
        raise type(error)(f"TidalPy: {where}: {detail}") from error
