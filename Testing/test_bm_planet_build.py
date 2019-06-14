import TidalPy
import numpy as np

planet = TidalPy.build_planet('charon_example', force_build=True)
host = TidalPy.build_planet('pluto_example', force_build=True)
star = TidalPy.build_planet('sol', force_build=True)

planet.paint(auto_show=True)

planet.mantle.temperature = 250.
planet.core.temperature = np.linspace(1100., 1700., 4)

orbit = TidalPy.Orbit(star, host, [planet])
orbit.set_eccentricity(2, .2)

print(planet.tidal_heating)
print('E', planet.core.tidal_heating)
print('eta', planet.core.viscosity)
print('mu', planet.core.shear_modulus)
print('k', planet.core.thermal_conductivity)
print('alpha', planet.core.thermal_expansion)
print('Rayleigh', planet.core.rayleigh)
print('nusselt', planet.core.nusselt)
print('blt', planet.core.blt/1e3)
print('heat_flux', planet.core.heat_flux)

planet.time = 0.

print(planet.layers_byname['core'].radiogenics.calculate())
print(planet.layers_byname['core'].calc_temperature_derivative())

print(planet.surface_temperature)