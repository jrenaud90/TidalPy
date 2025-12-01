import numpy as np
import TidalPy
TidalPy.config['stream_level'] = 'WARNING'
TidalPy.reinit()

from TidalPy.structures import build_world, build_from_world, Orbit

from performance_base import PerformanceTrackBase

star = build_world('55cnc')
io_base = build_world('io_simple')

class TideCalcPerformance(PerformanceTrackBase):

    def run_perform_build_calc_tides_one_layer(self):

        eccen_trunc = 6
        order_l = 2
        io = build_from_world(io_base, new_config={'tides': {'eccentricity_truncation_lvl': eccen_trunc,
                                                             'max_tidal_order_l': order_l,
                                                             'obliquity_tides_on': True,}})
        rheo = io.mantle.rheology.complex_compliance_model.model
        orbit = Orbit(star, tidal_host=star, tidal_bodies=io)
        io.eccentricity = 0.1
        io.obliquity = 0.1
        orbit.time = 0.
        set_temp_func = lambda: setattr(io.mantle, 'temperature', 1500.)

        self.record_performance(f'Calc Tides - One Layer - {rheo.title()} - Float', set_temp_func,
                                inputs=tuple(), repeats=10, number=500, note=f'Eccentricity Truc: {eccen_trunc}, '
                                                                              f'Tidal Order l: {order_l}.')

    def run_perform_build_calc_tides_one_layer_array(self):
        eccen_trunc = 6
        order_l = 2
        io = build_from_world(io_base, new_config={'tides': {'eccentricity_truncation_lvl': eccen_trunc,
                                                             'max_tidal_order_l'          : order_l,
                                                             'obliquity_tides_on'         : True, }})
        rheo = io.mantle.rheology.complex_compliance_model.model
        orbit = Orbit(star, tidal_host=star, tidal_bodies=io)
        io.eccentricity = 0.1
        io.obliquity = 0.1
        orbit.time = 0.
        temperature = np.linspace(1100., 2000., 10000)
        set_temp_func = lambda: setattr(io.mantle, 'temperature', temperature)

        self.record_performance(f'Calc Tides - One Layer - {rheo.title()} - Array', set_temp_func,
                                inputs=tuple(), repeats=10, number=300, array_N=len(temperature),
                                note=f'Eccentricity Truc: {eccen_trunc}, Tidal Order l: {order_l}.')


if __name__ == '__main__':
    TideCalcPerformance()