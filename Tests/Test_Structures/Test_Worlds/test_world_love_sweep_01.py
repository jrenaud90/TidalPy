"""Love-number sweeps over frequency (calc_love_numbers) and the quality factor and lag of the last solve."""
import math

import numpy as np
import pytest

FREQUENCIES = np.array([1.0e-5, 4.1e-5, 1.0e-4])


def test_sweep_matches_single_solves(io):
    sweep = io.calc_love_numbers(FREQUENCIES, degree_l=2)
    np.testing.assert_array_equal(sweep["frequency"], FREQUENCIES)
    assert sweep["success"].all()
    for frequency_i, frequency in enumerate(FREQUENCIES):
        single = io.solve_love_numbers(frequency=frequency, degree_l=2)
        assert sweep["k"][frequency_i] == single["love_number_k"]
        assert sweep["h"][frequency_i] == single["love_number_h"]
        assert sweep["l"][frequency_i] == single["love_number_l"]
        assert sweep["message"][frequency_i] == single["message"]


def test_sweep_keeps_the_input_shape(io):
    frequencies = FREQUENCIES.reshape(3, 1)
    sweep = io.calc_love_numbers(frequencies)
    assert sweep["k"].shape == (3, 1)
    assert sweep["success"].dtype == np.bool_
    scalar = io.calc_love_numbers(4.1e-5)
    assert scalar["k"].shape == ()


def test_failed_solves_are_nan_and_do_not_raise(io):
    sweep = io.calc_love_numbers(FREQUENCIES, max_num_steps=2)
    assert not sweep["success"].any()
    assert np.isnan(sweep["k"]).all()
    assert all(message for message in sweep["message"])


def test_sweep_refuses_its_own_arguments(io):
    with pytest.raises(TypeError, match="frequency"):
        io.calc_love_numbers(FREQUENCIES, frequency=1.0e-5)
    with pytest.raises(TypeError, match="raise_on_fail"):
        io.calc_love_numbers(FREQUENCIES, raise_on_fail=True)


def test_solve_love_numbers_can_raise_on_failure(io):
    from TidalPy.exceptions import SolutionFailedError
    with pytest.raises(RuntimeError, match="Love-number solve of world 'Io'"):
        io.solve_love_numbers(frequency=4.1e-5, max_num_steps=2, raise_on_fail=True)
    with pytest.raises(SolutionFailedError):
        io.solve_love_numbers(frequency=4.1e-5, max_num_steps=2, raise_on_fail=True)
    assert io.solve_love_numbers(frequency=4.1e-5, raise_on_fail=True)["success"]


@pytest.mark.parametrize("love_method", ["radial_solver", "homogeneous"])
def test_quality_factor_and_lag_of_k(io, love_method):
    io.solve_love_numbers(frequency=4.1e-5, love_method=love_method)
    love_k = io.love_number_k
    sign = -1.0 if love_k.real < 0.0 else 1.0
    assert io.love_q_k == pytest.approx(-sign * abs(love_k) / love_k.imag, rel=1.0e-14)
    assert io.love_lag_k == pytest.approx(math.atan2(-sign * love_k.imag, abs(love_k.real)), rel=1.0e-14)


def test_quality_factor_and_lag_are_nan_without_a_solve(io):
    unsolved = io.copy()
    assert math.isnan(unsolved.love_q_k)
    assert math.isnan(unsolved.love_lag_k)
