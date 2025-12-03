from performance_build_world import BuildWorldPerformance
from performance_complex_compliance_func import ComplexCompliancePerformance
from performance_eccentricity_func import EccentricityFuncPerformance
from performance_tides import TideCalcPerformance
from performance_quick_calcs import QuickCalcPerformance
from multilayer_radial_solver import TidalYPerformance
from multimode_solver import MultilayerModeNumbaPerformance, MultilayerModePerformance

if __name__ == '__main__':
    BuildWorldPerformance()
    ComplexCompliancePerformance()
    EccentricityFuncPerformance()
    TideCalcPerformance()
    QuickCalcPerformance()
    TidalYPerformance()
    MultilayerModeNumbaPerformance()
    MultilayerModePerformance()
