
import TidalPy

example_planet_config = TidalPy
example_planet = TidalPy.build_planet('earth_example', force_build=True, auto_save=False)

example_planet.upper_mantle.temperature = 1500.