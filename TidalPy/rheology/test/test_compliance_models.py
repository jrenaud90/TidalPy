""" Unit tests for rheology/compliance_models.py

Be sure to add a new test when you implement a new compliance model!
"""

import numpy as np

from ...types import float_eps
from ..compliance_models import off, fixed_q, maxwell, voigt, burgers, andrade, sundberg, andrade_freq, sundberg_freq

COMPLIANCE = 5.e10**-1
VISCOSITY = 1.e20
FREQUENCY = 2. * np.pi / (2. * 24. * 60. * 60.)

COMPLIANCE_ARR = np.asarray([COMPLIANCE])
VISCOSITY_ARR = np.asarray([VISCOSITY])
FREQUENCY_ARR = np.asarray([FREQUENCY])
NAN_ARR = np.asarray([np.nan])


class Test_ComplianceModels():

    @staticmethod
    def general_compliance_check(comp_func):

        # Basic NaN checks
        assert not np.isnan(np.real(comp_func(COMPLIANCE, VISCOSITY, FREQUENCY)))
        assert not np.isnan(np.imag(comp_func(COMPLIANCE, VISCOSITY, FREQUENCY)))
        assert not np.any(np.isnan(np.real(comp_func(COMPLIANCE_ARR, VISCOSITY_ARR, FREQUENCY_ARR))))
        assert not np.any(np.isnan(np.imag(comp_func(COMPLIANCE_ARR, VISCOSITY_ARR, FREQUENCY_ARR))))

        # Zero Frequency Check: Compliance functions should always return J + 0i for a zero frequency input.
        # TODO: Is that true for real(J^bar) in other rheologies?
        assert comp_func(COMPLIANCE, VISCOSITY, 0.0) == COMPLIANCE + 0.0j
        assert comp_func(COMPLIANCE_ARR, VISCOSITY_ARR, np.asarray([0.0])) == COMPLIANCE_ARR + \
               np.zeros(COMPLIANCE_ARR.shape, dtype=np.complex)
        assert comp_func(COMPLIANCE, VISCOSITY, float_eps) == COMPLIANCE + 0.0j
        assert comp_func(COMPLIANCE_ARR, VISCOSITY_ARR, np.asarray([float_eps])) == COMPLIANCE_ARR + \
               np.zeros(COMPLIANCE_ARR.shape, dtype=np.complex)
        assert comp_func(COMPLIANCE, VISCOSITY, -float_eps) == COMPLIANCE + 0.0j
        assert comp_func(COMPLIANCE_ARR, VISCOSITY_ARR, np.asarray([-float_eps])) == COMPLIANCE_ARR + \
               np.zeros(COMPLIANCE_ARR.shape, dtype=np.complex)

        # Check Bad Inputs
        assert np.isnan(np.real(comp_func(np.nan, VISCOSITY, FREQUENCY)))
        assert np.isnan(np.real(comp_func(np.nan, np.nan, FREQUENCY)))

        return True

    def test_off(self):

        # Check Floats
        float_result = off(COMPLIANCE, VISCOSITY, FREQUENCY)
        assert np.real(float_result) == COMPLIANCE
        assert np.imag(float_result) == 0.0
        assert float_result == COMPLIANCE + 0.0j

        # Check Arrays
        arr_result = off(COMPLIANCE_ARR, VISCOSITY_ARR, FREQUENCY_ARR)
        assert np.real(arr_result) == COMPLIANCE_ARR
        assert np.imag(arr_result) == np.asarray([0.])
        assert arr_result == COMPLIANCE_ARR + np.zeros(COMPLIANCE_ARR.shape, dtype=np.complex)

        # Check Bad Inputs
        assert np.imag(off(np.nan, VISCOSITY, FREQUENCY)) == 0.0
        assert np.imag(off(np.asarray([np.nan]), VISCOSITY_ARR, FREQUENCY_ARR)) == np.zeros_like(FREQUENCY_ARR)
        assert off(COMPLIANCE, np.nan, FREQUENCY) == COMPLIANCE + 0.0j
        assert off(COMPLIANCE_ARR, np.asarray([np.nan]), FREQUENCY_ARR) == COMPLIANCE_ARR + \
               np.zeros(FREQUENCY_ARR.shape, dtype=np.complex)

        assert self.general_compliance_check(off)
    #
    # def test_fixed_q(self):
    #
    #     # Check Floats
    #     float_result = fixed_q(COMPLIANCE_ARR, VISCOSITY_ARR, FREQUENCY_ARR, 10., 20.)
    #     assert np.testing.assert_approx_equal(np.real(float_result), np.asarray([-0.475]))
    #     assert np.testing.assert_approx_equal(np.imag(float_result), np.asarray([-112812500004.75003]))
    #
    #     assert self.general_compliance_check(fixed_q)