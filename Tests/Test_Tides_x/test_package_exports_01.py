"""The ``TidalPy.Tides_x`` package re-exports its public entry points; each is the object its subpackage holds."""
import importlib

import TidalPy.Tides_x as tides_x

_HOMES = {
    "TidalPy.Tides_x.classes": (
        "TideBase", "RheologyTide", "FixedQTide", "FixedLagTide", "CTLQTide", "make_tide", "collapse_global_tides"),
    "TidalPy.Tides_x.love": (
        "LoveNumbers", "apply_fixed_dt", "apply_fixed_q", "calc_effective_rigidity",
        "calc_homogeneous_love_numbers", "love_method_name"),
    "TidalPy.Tides_x.potential": ("ModeMap", "UniqueFrequencyMap", "tidal_potential_3d_modes", "global_potential"),
    "TidalPy.Tides_x.multilayer": (
        "angular_gram", "displacement_point", "strain_stress_heating_point", "volumetric_heating"),
    "TidalPy.Tides_x.eccentricity": ("eccentricity_func",),
    "TidalPy.Tides_x.obliquity": ("obliquity_func",),
}


def test_every_export_is_its_subpackage_object():
    for home, names in _HOMES.items():
        module = importlib.import_module(home)
        for name in names:
            assert name in tides_x.__all__, name
            assert getattr(tides_x, name) is getattr(module, name), name


def test_all_lists_exactly_the_exports():
    expected = {name for names in _HOMES.values() for name in names}
    assert set(tides_x.__all__) == expected
    assert len(tides_x.__all__) == len(expected)


def test_the_multilayer_subpackage_exports_the_kernels():
    from TidalPy.Tides_x import multilayer
    from TidalPy.Tides_x.multilayer import stress_strain
    for name in _HOMES["TidalPy.Tides_x.multilayer"]:
        assert getattr(multilayer, name) is getattr(stress_strain, name)
