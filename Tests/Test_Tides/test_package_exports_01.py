"""``TidalPy.Tides`` re-exports each public entry point as the object its subpackage holds."""
import importlib

import TidalPy.Tides as tides_module

_HOMES = {
    "TidalPy.Tides.classes": (
        "TideBase", "RheologyTide", "FixedQTide", "FixedLagTide", "CTLQTide", "make_tide", "collapse_global_tides"),
    "TidalPy.Tides.love": (
        "LoveNumbers", "apply_fixed_dt", "apply_fixed_q", "calc_effective_rigidity",
        "calc_homogeneous_love_numbers", "love_method_name"),
    "TidalPy.Tides.potential": ("ModeMap", "UniqueFrequencyMap", "tidal_potential_3d_modes", "global_potential"),
    "TidalPy.Tides.multilayer": (
        "angular_gram", "displacement_point", "strain_stress_heating_point", "volumetric_heating"),
    "TidalPy.Tides.eccentricity": ("eccentricity_func",),
    "TidalPy.Tides.obliquity": ("obliquity_func",),
}


def test_every_export_is_its_subpackage_object():
    """Each export is in __all__ and is the subpackage's object."""
    for home, names in _HOMES.items():
        module = importlib.import_module(home)
        for name in names:
            assert name in tides_module.__all__, name
            assert getattr(tides_module, name) is getattr(module, name), name


def test_all_lists_exactly_the_exports():
    """__all__ holds exactly the exports, without duplicates."""
    expected = {name for names in _HOMES.values() for name in names}
    assert set(tides_module.__all__) == expected
    assert len(tides_module.__all__) == len(expected)


def test_the_multilayer_subpackage_exports_the_kernels():
    """The multilayer subpackage re-exports the stress_strain kernels."""
    from TidalPy.Tides import multilayer
    from TidalPy.Tides.multilayer import stress_strain
    for name in _HOMES["TidalPy.Tides.multilayer"]:
        assert getattr(multilayer, name) is getattr(stress_strain, name)
