import TidalPy
import numpy as np

planet = TidalPy.build_planet('pluto_example', force_build=True)
host = TidalPy.build_planet('pluto_example', force_build=True)
star = TidalPy.build_planet('sol', force_build=True)

planet.paint(auto_show=True)

planet.core.config['surface_temperature'] = 275.15
planet.core.temperature = np.linspace(1100., 2000., 4)

orbit = TidalPy.Orbit(star, host, [planet])
orbit.set_semi_major_axis(2, TidalPy.conversions.Au2m(0.1))
orbit.set_eccentricity(2, .2)

print(planet.tidal_heating)

