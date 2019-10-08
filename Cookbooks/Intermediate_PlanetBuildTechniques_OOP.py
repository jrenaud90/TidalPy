from TidalPy.planets import build_planet, planet_iterator

# First lets build a reference object
trappist_1b = build_planet('trappist1b')
trappist_1b.paint()

new_trappist_1b = planet_iterator(trappist_1b, goal_radius=7.371e6, goal_mass=6.0737e24)
new_trappist_1b.paint()