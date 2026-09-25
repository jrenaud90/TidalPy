"""The default config's LLRI_and_SLRI isotope dataset (Castillo-Rogez et al. 2007).

The dataset once multiplied the paper's isotope concentrations by the isotopic abundances a second time, referenced
its formation-epoch abundances to 4600 Myr so that they applied today, and misread the 60Fe/56Fe ratio. Evaluated
today it returned 4.4e-12 W/kg from the long-extinct 26Al, and evaluated at formation it overflowed to infinity. The
dataset now reproduces the paper's Table 3 with time measured from formation.
"""
import math

import pytest
import toml

from TidalPy import defaultc
from TidalPy.radiogenics import radiogenic_models

SHORT_LIVED = ("Al26", "Fe60", "Mn53")


def _dataset(isotopes=None):
    known_isotope_data = toml.loads(defaultc.default_config_str)["physics"]["radiogenics"]["known_isotope_data"]
    dataset = known_isotope_data["LLRI_and_SLRI"]
    names = [name for name, value in dataset.items() if isinstance(value, dict)]
    if isotopes is not None:
        names = [name for name in names if name in isotopes]
    return (
        tuple(dataset[name]["iso_mass_fraction"] for name in names),
        tuple(dataset[name]["element_concentration"] for name in names),
        tuple(dataset[name]["half_life"] for name in names),
        tuple(dataset[name]["hpr"] for name in names),
        dataset["ref_time"])


def _specific_heating(time_myr, isotopes=None):
    """Heating per kilogram of rock [W kg-1]; the dataset's half lives and reference time are in Myr."""
    massfracs, concentrations, half_lives, hpr, ref_time = _dataset(isotopes)
    return float(radiogenic_models.isotope(time_myr, 1.0, massfracs, concentrations, half_lives, hpr, ref_time))


def test_reference_time_is_formation():
    assert _dataset()[4] == 0.0


def test_heating_at_formation():
    """Rock at formation heats at about 2.2e-7 W/kg, nearly all of it from 26Al."""
    heating = _specific_heating(0.0)
    assert math.isfinite(heating)
    assert heating == pytest.approx(2.2e-7, rel=0.05)
    assert _specific_heating(0.0, isotopes=("Al26",)) / heating > 0.9


def test_heating_today_is_long_lived():
    """At 4568 Myr the short-lived isotopes are extinct and the long-lived ones give a few 1e-12 W/kg."""
    heating = _specific_heating(4568.0)
    assert 4.0e-12 < heating < 6.0e-12
    assert _specific_heating(4568.0, isotopes=SHORT_LIVED) < 1.0e-20 * heating


@pytest.mark.parametrize("dataset_name", ("LLRI_and_SLRI", "llri_and_slri"))
def test_dataset_selected_by_name_in_a_world(dataset_name):
    """A layer's isotope dataset is found by name whatever its case (the mixed-case name once always raised)."""
    from TidalPy.structures import build_world, build_from_world

    io_base = build_world('io_simple')
    radiogenics = {'model': 'isotope', 'isotopes': dataset_name}
    radio_dict = {'layers': {'Mantle': {'is_tidal': False, 'radiogenics': radiogenics}}}
    io_llri = build_from_world(io_base, radio_dict)
    io_llri.set_state(time=0.)
    assert io_llri.Mantle.radiogenics.heating / io_llri.Mantle.mass == pytest.approx(2.2e-7, rel=0.05)


def test_unknown_dataset_name_raises():
    from TidalPy.exceptions import UnknownModelError
    from TidalPy.structures import build_world, build_from_world

    io_base = build_world('io_simple')
    radiogenics = {'model': 'isotope', 'isotopes': 'no_such_set'}
    radio_dict = {'layers': {'Mantle': {'is_tidal': False, 'radiogenics': radiogenics}}}
    with pytest.raises(UnknownModelError, match="no_such_set"):
        io_unknown = build_from_world(io_base, radio_dict)
        io_unknown.set_state(time=0.)
