"""A tide model built through its class takes the [tides] defaults for the lists it is not given, as make_tide does."""
import pytest

import TidalPy
from TidalPy.Tides.classes.tide import CTLQTide, FixedLagTide, FixedQTide, make_tide

FREQUENCY = 4.1e-5
GIVEN_K = [0.3]


@pytest.mark.parametrize("tide_class, model_name, keyword", [
    (FixedQTide, "cpl", "fixed_q"),
    (FixedLagTide, "ctl", "fixed_dt"),
    (CTLQTide, "ctl_q", "fixed_dt")])
def test_a_class_given_only_k_matches_make_tide(tide_class, model_name, keyword):
    by_class = tide_class(GIVEN_K)
    by_name = make_tide(model_name, {"fixed_k": GIVEN_K})
    expected = by_name.calc_neg_imk(2, FREQUENCY)
    assert expected > 0.0
    assert by_class.calc_neg_imk(2, FREQUENCY) == pytest.approx(expected, rel=1e-15)


def test_the_defaults_come_from_the_tides_config():
    tides = TidalPy.config["tides"]
    model = FixedQTide()
    assert model.get_fixed_k(2) == pytest.approx(tides["fixed_k"][0])
    assert model.get_fixed_q(2) == pytest.approx(tides["fixed_q"][0])
    assert FixedLagTide().get_fixed_dt(2) == pytest.approx(tides["fixed_dt_s"][0])


def test_given_lists_override_the_defaults():
    model = CTLQTide(fixed_k=[0.2], fixed_dt=[300.0], fixed_q=[40.0])
    assert model.get_fixed_k(2) == pytest.approx(0.2)
    assert model.get_fixed_dt(2) == pytest.approx(300.0)
    assert model.get_fixed_q(2) == pytest.approx(40.0)
