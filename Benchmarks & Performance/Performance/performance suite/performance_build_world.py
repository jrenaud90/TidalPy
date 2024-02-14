import TidalPy

TidalPy.config['stream_level'] = 'WARNING'
TidalPy.reinit()

from TidalPy.structures import build_world, build_from_world, Orbit


from performance_base import PerformanceTrackBase

class BuildWorldPerformance(PerformanceTrackBase):

    def run_perform_build_star(self):
        self.record_performance('Build World - Star', build_world,
                                inputs=('55cnc',), repeats=3, number=10)

    def run_perform_build_simple_layered(self):
        self.record_performance('Build World - Layered', build_world,
                                inputs=('io_simple',), repeats=3, number=10)


if __name__ == '__main__':
    performance_tracker = BuildWorldPerformance()