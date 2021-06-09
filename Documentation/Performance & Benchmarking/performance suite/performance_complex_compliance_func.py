import numpy as np

from performance_base import PerformanceTrackBase

import TidalPy
TidalPy.config['stream_level'] = 'WARNING'
TidalPy.reinit()

class ComplexCompliancePerformance(PerformanceTrackBase):

    def run_perform_complex_comp_maxwell_float(self):
        from TidalPy.rheology.complex_compliance.compliance_models import maxwell
        self.record_performance('Complex Compliance - Maxwell - Float', maxwell,
                                inputs=(1.1e-4, 2.1e-11, 3.3e9))

    def run_perform_complex_comp_andrade_float(self):
        from TidalPy.rheology.complex_compliance.compliance_models import andrade
        self.record_performance('Complex Compliance - Andrade - Float', andrade,
                                inputs=(1.1e-4, 2.1e-11, 3.3e9))

    def run_perform_complex_comp_sundberg_float(self):
        from TidalPy.rheology.complex_compliance.compliance_models import sundberg
        self.record_performance('Complex Compliance - Sundberg - Float', sundberg,
                                inputs=(1.1e-4, 2.1e-11, 3.3e9))

    def run_perform_complex_comp_maxwell_array(self):
        from TidalPy.rheology.complex_compliance.compliance_models import maxwell
        freq = np.linspace(1.1e-4, 5.1e-4, 10000)
        self.record_performance('Complex Compliance - Maxwell - Array', maxwell,
                                inputs=(freq, 2.1e-11, 3.3e9), array_N=len(freq))

    def run_perform_complex_comp_andrade_array(self):
        from TidalPy.rheology.complex_compliance.compliance_models import andrade
        freq = np.linspace(1.1e-4, 5.1e-4, 10000)
        self.record_performance('Complex Compliance - Andrade - Array', andrade,
                                inputs=(freq, 2.1e-11, 3.3e9), array_N=len(freq))

    def run_perform_complex_comp_sundberg_array(self):
        from TidalPy.rheology.complex_compliance.compliance_models import sundberg
        freq = np.linspace(1.1e-4, 5.1e-4, 10000)
        self.record_performance('Complex Compliance - Sundberg - Array', sundberg,
                                inputs=(freq, 2.1e-11, 3.3e9), array_N=len(freq))


if __name__ == '__main__':
    performance_tracker = ComplexCompliancePerformance()