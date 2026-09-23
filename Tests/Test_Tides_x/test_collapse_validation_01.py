"""The standalone ``collapse_global_tides`` checks its input the way ``make_tide`` and a world do.

An absent config takes the ``[tides]`` defaults instead of giving no heating, a misspelled key raises, a whole-valued
float truncation is the level it names, and an out-of-order degree range or an eccentricity outside [0, 1) raises.
"""
import pytest

from TidalPy.constants import G
from TidalPy.Tides_x.classes.collapse import collapse_global_tides

_IO = dict(planet_radius=1.8215e6, orbital_frequency=4.11e-5, spin_frequency=4.11e-5, eccentricity=0.0041,
           obliquity=0.001, semi_major_axis=4.217e8, host_mass=1.898e27, G_to_use=G)


def test_an_absent_config_takes_the_defaults():
    result = collapse_global_tides(**_IO, tide_model="cpl")
    assert result["tidal_heating"] > 0.0


def test_a_misspelled_key_raises():
    with pytest.raises(ValueError):
        collapse_global_tides(**_IO, tide_model="ctl", tide_config={"fixd_dt_s": [100.0]})


def test_a_whole_valued_float_truncation_is_its_level():
    config = {"fixed_k": [0.3], "fixed_q": [100.0]}
    as_int = collapse_global_tides(**_IO, tide_model="cpl", tide_config=config, obliquity_truncation=2)
    as_float = collapse_global_tides(**_IO, tide_model="cpl", tide_config=config, obliquity_truncation=2.0)
    assert as_float["tidal_heating"] == as_int["tidal_heating"]
    with pytest.raises(ValueError):
        collapse_global_tides(**_IO, tide_model="cpl", tide_config=config, obliquity_truncation=2.5)


@pytest.mark.parametrize("override", [dict(eccentricity=1.5), dict(eccentricity=-0.1), dict(semi_major_axis=0.0)])
def test_an_impossible_orbit_raises(override):
    with pytest.raises(ValueError):
        collapse_global_tides(**dict(_IO, **override), tide_model="cpl")


def test_an_out_of_order_degree_range_raises():
    with pytest.raises(ValueError):
        collapse_global_tides(**_IO, tide_model="cpl", min_degree_l=3, max_degree_l=2)
