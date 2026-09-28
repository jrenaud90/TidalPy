"""Input checks of the standalone ``collapse_global_tides``, matching ``make_tide`` and the world."""
import pytest

from TidalPy.constants import G
from TidalPy.Tides.classes.collapse import collapse_global_tides

_IO = dict(planet_radius=1.8215e6, orbital_frequency=4.11e-5, spin_frequency=4.11e-5, eccentricity=0.0041,
           obliquity=0.001, semi_major_axis=4.217e8, host_mass=1.898e27, G_to_use=G)


def test_an_absent_config_takes_the_defaults():
    """No config uses the [tides] defaults rather than giving no heating."""
    result = collapse_global_tides(**_IO, tide_model="cpl")
    assert result["tidal_heating"] > 0.0


def test_a_whole_valued_float_truncation_is_its_level():
    """A truncation of 2.0 is level 2; 2.5 raises."""
    config = {"fixed_k": [0.3], "fixed_q": [100.0]}
    as_int = collapse_global_tides(**_IO, tide_model="cpl", tide_config=config, obliquity_truncation=2)
    as_float = collapse_global_tides(**_IO, tide_model="cpl", tide_config=config, obliquity_truncation=2.0)
    assert as_float["tidal_heating"] == as_int["tidal_heating"]
    with pytest.raises(TypeError):
        collapse_global_tides(**_IO, tide_model="cpl", tide_config=config, obliquity_truncation=2.5)


@pytest.mark.parametrize("override", [
    pytest.param(dict(eccentricity=1.5), id="eccentricity_above_one"),
    pytest.param(dict(eccentricity=-0.1), id="negative_eccentricity"),
    pytest.param(dict(semi_major_axis=0.0), id="zero_semi_major_axis"),
    pytest.param(dict(tide_model="ctl", tide_config={"fixd_dt_s": [100.0]}), id="misspelled_key"),
    pytest.param(dict(min_degree_l=3, max_degree_l=2), id="out_of_order_degrees"),
])
def test_invalid_input_raises(override):
    """An impossible orbit, a misspelled config key, or an out-of-order degree range raises."""
    with pytest.raises(ValueError):
        collapse_global_tides(**{**_IO, "tide_model": "cpl", **override})
