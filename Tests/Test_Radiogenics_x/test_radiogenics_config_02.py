"""A layer's own isotope arrays win over a dataset name, and a given reference time applies to a dataset."""
import pytest

from TidalPy.radiogenics_x.radiogenics import make_radiogenics

MYR = 3.15576e13  # [s]
EXPLICIT = {
    "heat_production_w_kg": [1.0e-5],
    "half_lives_s": [1.0e30],
    "mass_fracs": [1.0],
    "concentrations": [1.0],
}


def test_explicit_arrays_win_over_a_dataset_name():
    alone = make_radiogenics("isotope", dict(EXPLICIT))
    with_name = make_radiogenics("isotope", {**EXPLICIT, "isotopes": "modern_day_chondritic"})
    assert with_name.calc_heating(0.0, 1.0) == pytest.approx(alone.calc_heating(0.0, 1.0), rel=1e-12)
    chondritic = make_radiogenics("isotope", {"isotopes": "modern_day_chondritic"})
    assert with_name.calc_heating(0.0, 1.0) != pytest.approx(chondritic.calc_heating(0.0, 1.0), rel=1e-3)


def test_reference_time_applies_to_a_dataset():
    default = make_radiogenics("isotope", {"isotopes": "llri_and_slri"})
    shifted = make_radiogenics("isotope", {"isotopes": "llri_and_slri", "ref_time_s": 100.0 * MYR})
    # The same abundances placed 100 Myr later give, at 100 Myr, the heating the dataset gives at its own epoch.
    assert shifted.calc_heating(100.0 * MYR, 1.0) == pytest.approx(default.calc_heating(0.0, 1.0), rel=1e-9)
