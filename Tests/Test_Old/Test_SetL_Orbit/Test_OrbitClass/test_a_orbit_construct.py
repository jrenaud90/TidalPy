import copy

import TidalPy


from TidalPy.structures import build_world, scale_from_world
from TidalPy.structures.orbit import OrbitBase

big_io_config = {
    "name"           : "BigIo",
    "type"           : "layered",
    "radius"         : 1821.49e3 * 10.,
    "orbital_period" : 1.769,
    "eccentricity"   : 0.0041,
    "spin_period"    : 1.769,
    "albedo"         : 0.63,
    "force_spin_sync": True,
    "layers"         : {
        "Core"  : {
            "type"    : "iron",
            "is_tidal": False,
            "radius"  : 810.0e3,
            "density" : 5200.
            },
        "Mantle": {
            "type"               : "rock",
            "is_tidal"           : True,
            "radius"             : 1821.49e3 * 10.,
            "surface_temperature": 100.0,
            "density"            : 3200.
            }
        }
    }


def test_basic_orbit_construction():
    # Construct orbit with no members
    orbit = OrbitBase()

    assert type(orbit) == OrbitBase


def test_orbit_construction_with_star():
    """ This will test building an orbit with just a star. """

    star = build_world('55Cnc')
    star_2 = copy.deepcopy(star)

    # Simple check to see if star is loaded in correctly
    orbit_1 = OrbitBase(star=star)
    assert orbit_1.star is star
    assert orbit_1.all_objects[0] is star
    assert star.orbit is orbit_1

    # Check if star loaded as tidal host works as expected
    orbit_2 = OrbitBase(star=star_2, star_host=True)
    assert orbit_2.star is star_2
    assert orbit_2.all_objects[0] is star_2
    assert star.orbit is orbit_1
    assert orbit_2.tidal_host is star_2
    assert orbit_2.tidal_objects[0] is star_2


def test_orbit_construction_with_star_and_host():
    """ This will test building an orbit with just a star and tidal host. """

    star = build_world('55Cnc')
    host = build_world('io_simple')

    # Simple check to see if star and host are loaded in correctly
    orbit = OrbitBase(star=star, tidal_host=host)
    assert orbit.star is star
    assert orbit.tidal_host is host
    assert orbit.all_objects[0] is star
    assert orbit.all_objects[1] is host
    assert star.orbit is orbit
    assert host.orbit is orbit

    # Tidal host should be the first (and only) tidal body
    assert len(orbit.tidal_objects) == 1
    assert orbit.tidal_objects[0] is host


def test_orbit_construction_with_star_and_host_and_tidalbody():
    """ This will test building an orbit with one tidal world. """

    star = build_world('55Cnc')
    big_io = build_world('BigIo', world_config=big_io_config)
    io = build_world('io_simple')

    star_2 = copy.deepcopy(star)
    big_io_2 = copy.deepcopy(big_io)
    io_2 = copy.deepcopy(io)

    # Simple check to see if star, host, and tidal body are loaded in correctly.
    orbit = OrbitBase(star=star, tidal_host=big_io, tidal_bodies=io)
    assert orbit.star is star
    assert orbit.tidal_host is big_io
    assert orbit.all_objects[0] is star
    assert orbit.all_objects[1] is big_io
    assert orbit.all_objects[2] is io
    assert star.orbit is orbit
    assert big_io.orbit is orbit
    assert io.orbit is orbit

    # Tidal host should be the first tidal body followed by other tidal bodies
    assert len(orbit.tidal_objects) == 2
    assert orbit.tidal_objects[0] is big_io
    assert orbit.tidal_objects[1] is io

    # Check that signatures are working as intended. (the tidal host should return the tidal body's orbit index)
    assert orbit.world_signature_to_index(big_io) == 1
    assert orbit.world_signature_to_index(big_io.name) == 1
    assert orbit.world_signature_to_index(0) == 1
    assert orbit.world_signature_to_index(io) == 1
    assert orbit.world_signature_to_index(io.name) == 1
    assert orbit.world_signature_to_index(1) == 1

    # Check that everything works as expected if one tidal body is passed as a list.
    orbit_2 = OrbitBase(star=star_2, tidal_host=big_io_2, tidal_bodies=[io_2])
    assert orbit_2.star is star_2
    assert orbit_2.tidal_host is big_io_2
    assert orbit_2.all_objects[0] is star_2
    assert orbit_2.all_objects[1] is big_io_2
    assert orbit_2.all_objects[2] is io_2
    assert star_2.orbit is orbit_2
    assert big_io_2.orbit is orbit_2
    assert io_2.orbit is orbit_2

    # Tidal host should be the first tidal body followed by other tidal bodies
    assert len(orbit_2.tidal_objects) == 2
    assert orbit_2.tidal_objects[0] is big_io_2
    assert orbit_2.tidal_objects[1] is io_2


def test_orbit_construction_with_no_star_and_host_and_tidalbody():
    """ This will test building an orbit with one tidal world, but no star. """

    big_io = build_world('BigIo', world_config=big_io_config)
    io = build_world('io_simple')

    # Simple check to see if star, host, and tidal body are loaded in correctly.
    orbit = OrbitBase(tidal_host=big_io, tidal_bodies=io)
    assert orbit.tidal_host is big_io
    assert orbit.all_objects[0] is big_io
    assert orbit.all_objects[1] is io
    assert big_io.orbit is orbit
    assert io.orbit is orbit

    # Tidal host should be the first tidal body followed by other tidal bodies
    assert len(orbit.tidal_objects) == 2
    assert orbit.tidal_objects[0] is big_io
    assert orbit.tidal_objects[1] is io

    # Check that signatures are working as intended. (the tidal host should return the tidal body's orbit index)
    assert orbit.world_signature_to_index(big_io) == 1
    assert orbit.world_signature_to_index(big_io.name) == 1
    assert orbit.world_signature_to_index(0) == 1
    assert orbit.world_signature_to_index(io) == 1
    assert orbit.world_signature_to_index(io.name) == 1
    assert orbit.world_signature_to_index(1) == 1


def test_orbit_construction_with_star_and_host_and_multi_tidalbodies():
    """ This will test building an orbit with multiple tidal world_types. """

    star = build_world('55Cnc')
    big_io = build_world('BigIo', world_config=big_io_config)
    io = build_world('io_simple')
    # Build some additional world_types to use as tidal bodies
    small_io = scale_from_world(io, new_name='small_io', radius_scale=0.5)
    tiny_io = scale_from_world(io, new_name='tiny_io', radius_scale=0.1)

    # Construct orbit using all of the tidal bodies
    orbit = OrbitBase(star=star, tidal_host=big_io, tidal_bodies=[io, small_io, tiny_io])
    assert orbit.star is star
    assert orbit.tidal_host is big_io
    assert orbit.all_objects[0] is star
    assert orbit.all_objects[1] is big_io
    assert orbit.all_objects[2] is io
    assert orbit.all_objects[3] is small_io
    assert orbit.all_objects[4] is tiny_io
    assert star.orbit is orbit
    assert big_io.orbit is orbit
    assert io.orbit is orbit
    assert small_io.orbit is orbit
    assert tiny_io.orbit is orbit

    # Tidal host should be the first tidal body followed by other tidal bodies
    assert len(orbit.tidal_objects) == 4
    assert orbit.tidal_objects[0] is big_io
    assert orbit.tidal_objects[1] is io
    assert orbit.tidal_objects[2] is small_io
    assert orbit.tidal_objects[3] is tiny_io

    # Check that signatures are working as intended.
    assert orbit.world_signature_to_index(big_io) == 1
    assert orbit.world_signature_to_index(big_io.name) == 1
    assert orbit.world_signature_to_index(0) == 1
    assert orbit.world_signature_to_index(io) == 1
    assert orbit.world_signature_to_index(io.name) == 1
    assert orbit.world_signature_to_index(1) == 1
    assert orbit.world_signature_to_index(small_io) == 2
    assert orbit.world_signature_to_index(small_io.name) == 2
    assert orbit.world_signature_to_index(2) == 2
    assert orbit.world_signature_to_index(tiny_io) == 3
    assert orbit.world_signature_to_index(tiny_io.name) == 3
    assert orbit.world_signature_to_index(3) == 3
