from TidalPy import build_world, build_from_world
from TidalPy.orbit import PhysicsOrbit

star = build_world('55cnc')
world = build_world('earth_simple')

def test_apply_physics_orbit():
    """ This will test applying a Physical orbit to the system """

    orbit = PhysicsOrbit(star, tidal_host=star, tidal_bodies=world)

    assert orbit.tidal_host is star
    assert orbit.star is star
    assert orbit.tidal_objects[0] is star
    assert orbit.tidal_objects[1] is world
