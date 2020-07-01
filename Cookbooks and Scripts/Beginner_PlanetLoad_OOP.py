import copy
from TidalPy.planets import build_planet

# Build the Earth
earth = build_planet('earth')
earth.paint()

# We are going to need a second copy for later
earth_2 = copy.deepcopy(earth)

# Setup the orbit
#  In this example we will consider two orbits: Heliocentric and Geocentric
from TidalPy.orbit import Orbit
sun = build_planet('sol')
# The Sun is entered into the orbit twice because, for this orbit, we are calculating heliocentric tides
#  instead of lunar tides.
helio_orbit = Orbit(sun, sun, [earth])

# For the geocentric orbit we are actually going to set it up so that the moon is the 'host'.
moon = build_planet('luna')
moon.paint()
geo_orbit = Orbit(sun, moon, [earth_2], use_host_orbit=True)

# Earth's configuration file restricts tidal heating to the upper mantle.
earth.upper_mantle.temperature = 1600.
earth_2.upper_mantle.temperature = 1600.
# In reality, for the earth, tides will dominate in the liquid ocean. But, TidalPy is not designed to calculate
#  tidal dissipation in liquids.

# Now we have all the information to calculate tidal heating
from TidalPy.tools.conversions import m2Au
print(f'Earth-Sun eccentricity = {earth.eccentricity}')
print(f'Earth-Sun semi-major axis = {m2Au(earth.semi_major_axis)} Au')
print(f'Solar tidal heating = {(earth.tidal_heating/1e6)} MW')

print(f'Earth-Moon eccentricity = {earth_2.eccentricity}')
print(f'Earth-Moon semi-major axis = {m2Au(earth_2.semi_major_axis)} Au')
print(f'Lunar tidal heating = {(earth_2.tidal_heating/1e6)} MW')

