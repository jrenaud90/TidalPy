"""ModelFamily: the Python-level registry of one spec-driven physics family.

Every family built on parameter specs (viscosity, rheology, equations of state, ...) exposes the same module-level
functions: its model names, a name's canonical form, the config keys a model reads, an alias-aware name comparison
(MatPack's merge uses it), and a ``make_*`` factory. A ModelFamily holds the family's wrapper classes and provides those
functions once, so a family module only names its classes and wraps these in documented functions.
"""

# Family name (as the C++ models report it) -> its ModelFamily, filled as each family module is imported.
_FAMILIES = {}


def get_family(family: str) -> "ModelFamily":
    """A registered family by its name (``"viscosity"``, ``"equation of state"``, ...).

    Raises
    ------
    KeyError
        No family of that name has been registered (its module was not imported).
    """
    return _FAMILIES[family]


def model_class(family: str, model_name: str):
    """The Python class of a model, by its family and canonical model name.

    Raises
    ------
    TypeError
        No Python class is registered for that family and model.
    """
    family_record = _FAMILIES.get(family)
    if family_record is None or model_name not in family_record.classes:
        raise TypeError(f"TidalPy: no Python class is registered for the {family} model '{model_name}'.")
    return family_record.classes[model_name]


class ModelFamily:
    """The wrapper classes of one physics family and the lookups every family shares.

    Parameters
    ----------
    family : str
        Family name used in messages, for example ``"viscosity"``.
    classes : sequence of type
        The concrete wrapper classes; each has a ``MODEL_NAME`` attribute, its canonical model name.
    canonical_name : callable
        ``canonical_name(name) -> str``, the C++ registry's alias-aware, case-insensitive name lookup; raises
        ``ValueError`` for an unknown name.
    """

    def __init__(self, family, classes, canonical_name):
        self.family = family
        self.classes = {cls.MODEL_NAME: cls for cls in classes}
        self.canonical_name = canonical_name
        # Each model's config keys, read from a default instance's parameter table.
        # A composite class lists the keys of its sub-model tables in EXTRA_CONFIG_KEYS.
        self.model_config_keys = {
            name: frozenset(entry["key"] for entry in cls().get_parameter_info())
                  | frozenset(getattr(cls, "EXTRA_CONFIG_KEYS", ()))
            for name, cls in self.classes.items()}
        self.config_keys = frozenset().union(*self.model_config_keys.values())
        _FAMILIES[family] = self

    def model_names(self) -> tuple:
        """The canonical model names, in registry order."""
        return tuple(self.classes)

    def config_keys_of(self, model_name: str) -> frozenset:
        """The config keys one model reads, by any of its names."""
        return self.model_config_keys[self.canonical_name(model_name)]

    def same_model(self, table_name: str, model_name: str) -> bool:
        """Whether two names (aliases included) resolve to the same model."""
        return self.canonical_name(table_name) == self.canonical_name(model_name)

    def make(self, model_name: str, config=None):
        """Build a model by name from a config dict keyed by config key; absent keys (all of them for ``None``) take
        the model's defaults."""
        return self.classes[self.canonical_name(model_name)](config=config)
