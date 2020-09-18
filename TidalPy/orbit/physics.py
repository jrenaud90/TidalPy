from typing import Union, List

from .base import OrbitBase
from ..exceptions import InitiatedPropertyChangeError
from ..structures.worlds import AllWorldType, StarWorld, all_tidal_world_types

WorldSignatureType = Union[str, int, AllWorldType]


class PhysicsOrbit(OrbitBase):

    """ PhysicsOrbit class
    Contains attributes and methods to track the orbit of multiple TidalPy worlds. Also contains attributes and methods
        used to calculate various physics for the worlds stored within. This includes: tidal heating, tidal evolution,
        insolation heating, etc.

    Orbits allow TidalPy worlds to communicate with one another and for tides to be calculated.

    Assumptions
    -----------
    .. All TidalPy orbits currently assume no interaction between tidal bodies or the orbit's star. The only interaction
        that is permitted is between a single tidal body and the orbit's tidal host (which could be the star).

    See Also
    --------
    TidalPy.orbit.OrbitBase
    """

    def __init__(self, star: StarWorld = None, tidal_host: AllWorldType = None,
                 tidal_bodies: Union[AllWorldType, List[AllWorldType]] = None, star_host: bool = False,
                 host_tide_raiser: AllWorldType = None, initialize: bool = True):
        """ PhysicsOrbit constructor

        Orbits allow TidalPy worlds to communicate with one another and for tides to be calculated.

        Notes
        -----
        .. Orbit instances can be made without any TidalPy worlds initially. Worlds can be later added with the
            `add_tidal_world` and `add_star` methods.

        Parameters
        ----------
        star : StarWorld
            Initialized TidalPy StarWorld that will become this orbit's star.
        tidal_host : AllWorldType
            Initialized TidalPy world that will become this orbit's tidal host.
        tidal_bodies : Union[AllWorldType, List[AllWorldType]]
            Initialized TidalPy world(s) that will become connected via this orbit to the tidal host and star.
            Multiple tidal bodies can be entered at once if provided to the constructor as a list.
        star_host : bool = False
            If `True`, then the star will act as the tidal host in addition to being the orbit's star.
        host_tide_raiser : AllWorldType = None
            The tidal host experiences tides from all worlds. This pointer is used to set which tidal body is currently
                being used as the host body's tide raiser.
        initialize : bool = True
            If `True`, then the constructor will make the first call to the orbit's `reinit()` method.

        See Also
        --------
        TidalPy.orbit.OrbitBase
        """

        super().__init__(star, tidal_host, tidal_bodies,
                         star_host=star_host, host_tide_raiser=host_tide_raiser, initialize=False)

        # Initialized properties
        self._tidally_active_worlds = list()

        # State properties


        if initialize:
            # Make first call to the classes reinit method.
            self.reinit(initial_init=True, reinit_worlds=False, run_update=True)

    def add_tidal_world(self, tidal_world: AllWorldType, is_tidal_host: bool = False, run_update: bool = True):
        """ Add a new tidal world to the orbit, in order from closest to host to farthest away.

        Parameters
        ----------
        tidal_world : AllWorldType
            TidalPy world instance to be added to orbit.
        is_tidal_host : bool = False
            If true, then additional checks will be done to ensure proper functionality.
        run_update : bool = True
            If `True`, the orbit's `update_orbit` method will be called after the world has been added.
        """

        # Add the world via the parent class method.
        super().add_tidal_world(tidal_world, is_tidal_host=is_tidal_host, run_update=run_update)

        # If the world is tidal then add it to the tidally active list
        if isinstance(tidal_world, all_tidal_world_types):
            # Check if the world can even be tidally active or not.
            if tidal_world.tides_on:
                # Check if tides are turned on or not inside world.
                if tidal_world not in self.tidally_active_worlds:
                    self._tidally_active_worlds.append(tidal_world)


    def update_orbit(self):
        """ Method that updates anything that depends upon the state orbital properties of the tidal worlds.

        Should be called whenever there is a change to one or more orbital properties.
        """

        super().update_orbit()

        for tidally_active_world in self.tidally_active_worlds:
            tidally_active_world.update_orbit_spin(called_from_orbit=True)


    # # Initialized properties
    @property
    def tidally_active_worlds(self):
        """ List of worlds stored in orbit that are tidally active. """
        return self._tidally_active_worlds

    @tidally_active_worlds.setter
    def tidally_active_worlds(self, value):
        raise InitiatedPropertyChangeError