from TidalPy.planets import build_planet, planet_iterator
from TidalPy.constants import M_earth

# First lets build a reference object
trappist_1b = build_planet('trappist1b')
trappist_1b.paint()

trappist_1c = build_planet('trappist1c')
trappist_1c.paint()

trappist_1d = build_planet('trappist1d')
trappist_1d.paint()

trappist_1e = build_planet('trappist1e')
trappist_1e.paint()

trappist_1f = build_planet('trappist1f')
trappist_1f.paint()

trappist_1g = build_planet('trappist1g')
trappist_1g.paint()

print('TRAP-1b:', trappist_1b.mass / (1.017 * M_earth))
print('TRAP-1c:', trappist_1c.mass / (1.156 * M_earth))
print('TRAP-1d:', trappist_1d.mass / (0.297 * M_earth))
print('TRAP-1e:', trappist_1e.mass / (0.772 * M_earth))
print('TRAP-1f:', trappist_1f.mass / (0.934 * M_earth))
print('TRAP-1g:', trappist_1g.mass / (1.148 * M_earth))


#new_trappist_1b = planet_iterator(trappist_1b, goal_radius=7.371e6, goal_mass=6.0737e24)
#new_trappist_1b.paint()