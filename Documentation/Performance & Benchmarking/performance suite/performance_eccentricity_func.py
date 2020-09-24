import numpy as np

from performance_base import PerformanceTrackBase

import TidalPy
TidalPy.config['stream_level'] = 'WARNING'
TidalPy.reinit()

class EccentricityFuncPerformance(PerformanceTrackBase):

    def run_perform_eccentricity_l2e2_float(self):
        from TidalPy.tides.eccentricity_funcs import eccentricity_funcs_l2_trunc2
        self.record_performance('Eccentricity Func (l=2, e^2) - Float', eccentricity_funcs_l2_trunc2,
                                inputs=(0.2,))

    def run_perform_eccentricity_l2e10_float(self):
        from TidalPy.tides.eccentricity_funcs import eccentricity_funcs_l2_trunc10
        self.record_performance('Eccentricity Func (l=2, e^10) - Float', eccentricity_funcs_l2_trunc10,
                                inputs=(0.2,))

    def run_perform_eccentricity_l2e20_float(self):
        from TidalPy.tides.eccentricity_funcs import eccentricity_funcs_l2_trunc20
        self.record_performance('Eccentricity Func (l=2, e^20) - Float', eccentricity_funcs_l2_trunc20,
                                inputs=(0.2,))

    def run_perform_eccentricity_l2e2_array(self):
        from TidalPy.tides.eccentricity_funcs import eccentricity_funcs_l2_trunc2
        eccentricity = np.linspace(0.1, 0.5, 10000)
        self.record_performance('Eccentricity Func (l=2, e^2) - Array', eccentricity_funcs_l2_trunc2,
                                inputs=(eccentricity,), array_N=len(eccentricity))

    def run_perform_eccentricity_l2e10_array(self):
        from TidalPy.tides.eccentricity_funcs import eccentricity_funcs_l2_trunc10
        eccentricity = np.linspace(0.1, 0.5, 10000)
        self.record_performance('Eccentricity Func (l=2, e^10) - Array', eccentricity_funcs_l2_trunc10,
                                inputs=(eccentricity,), array_N=len(eccentricity))

    def run_perform_eccentricity_l2e20_array(self):
        from TidalPy.tides.eccentricity_funcs import eccentricity_funcs_l2_trunc20
        eccentricity = np.linspace(0.1, 0.5, 10000)
        self.record_performance('Eccentricity Func (l=2, e^20) - Array', eccentricity_funcs_l2_trunc20,
                                inputs=(eccentricity,), array_N=len(eccentricity))

    def run_perform_eccentricity_l4e2_float(self):
        from TidalPy.tides.eccentricity_funcs import eccentricity_funcs_l2_trunc2
        self.record_performance('Eccentricity Func (l=4, e^2) - Float', eccentricity_funcs_l2_trunc2,
                                inputs=(0.2,))

    def run_perform_eccentricity_l4e10_float(self):
        from TidalPy.tides.eccentricity_funcs import eccentricity_funcs_l2_trunc10
        self.record_performance('Eccentricity Func (l=4, e^10) - Float', eccentricity_funcs_l2_trunc10,
                                inputs=(0.2,))

    def run_perform_eccentricity_l4e20_float(self):
        from TidalPy.tides.eccentricity_funcs import eccentricity_funcs_l2_trunc20
        self.record_performance('Eccentricity Func (l=4, e^20) - Float', eccentricity_funcs_l2_trunc20,
                                inputs=(0.2,))

    def run_perform_eccentricity_l4e2_array(self):
        from TidalPy.tides.eccentricity_funcs import eccentricity_funcs_l2_trunc2
        eccentricity = np.linspace(0.1, 0.5, 10000)
        self.record_performance('Eccentricity Func (l=4, e^2) - Array', eccentricity_funcs_l2_trunc2,
                                inputs=(eccentricity,), array_N=len(eccentricity))

    def run_perform_eccentricity_l4e10_array(self):
        from TidalPy.tides.eccentricity_funcs import eccentricity_funcs_l2_trunc10
        eccentricity = np.linspace(0.1, 0.5, 10000)
        self.record_performance('Eccentricity Func (l=4, e^10) - Array', eccentricity_funcs_l2_trunc10,
                                inputs=(eccentricity,), array_N=len(eccentricity))

    def run_perform_eccentricity_l4e20_array(self):
        from TidalPy.tides.eccentricity_funcs import eccentricity_funcs_l2_trunc20
        eccentricity = np.linspace(0.1, 0.5, 10000)
        self.record_performance('Eccentricity Func (l=4, e^20) - Array', eccentricity_funcs_l2_trunc20,
                                inputs=(eccentricity,), array_N=len(eccentricity))

if __name__ == '__main__':
    performance_tracker = EccentricityFuncPerformance()