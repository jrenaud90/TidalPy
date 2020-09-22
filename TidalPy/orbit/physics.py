from typing import Union, List

import numpy as np

from .base import OrbitBase, WorldSignatureType
from .. import log
from ..exceptions import InitiatedPropertyChangeError
from ..structures.world_types import AllWorldType, StarWorld, all_tidal_world_types, TidalWorldType
from ..utilities.types import FloatArray


class PhysicsOrbit(OrbitBase):

    """ PhysicsOrbit class
    Contains attributes and methods to track the orbit of multiple TidalPy world_types. Also contains attributes and methods
        used to calculate various physics for the world_types stored within. This includes: tidal heating, tidal evolution,
        insolation heating, etc.

    Orbits allow TidalPy world_types to communicate with one another and for tides to be calculated.

    Assumptions
    -----------
    .. All TidalPy orbits currently assume no interaction between tidal bodies or the orbit's star. The only interaction
        that is permitted is between a single tidal body and the orbit's tidal host (which could be the star).

    See Also
    --------
    TidalPy.orbit.OrbitBase
    """

    class_name = 'physics'

    def __init__(self, star: StarWorld = None, tidal_host: AllWorldType = None,
                 tidal_bodies: Union[AllWorldType, List[AllWorldType]] = None, star_host: bool = False,
                 host_tide_raiser: AllWorldType = None, initialize: bool = True):
        """ PhysicsOrbit constructor

        Orbits allow TidalPy world_types to communicate with one another and for tides to be calculated.

        Notes
        -----
        .. Orbit instances can be made without any TidalPy world_types initially. Worlds can be later added with the
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
            The tidal host experiences tides from all world_types. This pointer is used to set which tidal body is currently
                being used as the host body's tide raiser.
        initialize : bool = True
            If `True`, then the constructor will make the first call to the orbit's `reinit()` method.

        See Also
        --------
        TidalPy.orbit.OrbitBase
        """

        # State properties
        self._tidally_active_worlds = list()

        super().__init__(star, tidal_host, tidal_bodies,
                         star_host=star_host, host_tide_raiser=host_tide_raiser, initialize=False)

        if initialize:
            # Make first call to the methods reinit method.
            self.reinit(initial_init=True, reinit_worlds=False, run_update=True)

    def orbit_changed(self, specific_world: WorldSignatureType = None, orbital_freq_changed: bool = False,
                      eccentricity_changed: bool = False):
        """ The orbit of a specific world has changed. Make any necessary updates.

        Parameters
        ----------
        specific_world : WorldSignatureType
            The signature of the world who's spin and/or orbit changed.
        orbital_freq_changed : bool = False
            If `True`, then the world's orbital frequency changed.
        eccentricity_changed : bool = False
            If `True`, then the world_types' eccentricity changed.
        """

        super().orbit_changed(specific_world, orbital_freq_changed=orbital_freq_changed,
                              eccentricity_changed=eccentricity_changed)

        # Tell the world that its orbit has changed.
        world_index = self.world_signature_to_index(specific_world)
        world_instance = self.tidal_objects[world_index]

        # If world_instance is the host's tide raiser, then any orbital changes done to it will also effect the host.
        if world_instance is self.host_tide_raiser:
            # Update the host world first.
            self.tidal_host.orbit_spin_changed(orbital_freq_changed=orbital_freq_changed, spin_freq_changed=False,
                                               eccentricity_changed=eccentricity_changed, obliquity_changed=False,
                                               call_orbit_dissipation=False)
        # Now update the target world
        world_instance.orbit_spin_changed(orbital_freq_changed=orbital_freq_changed, spin_freq_changed=False,
                                          eccentricity_changed=eccentricity_changed, obliquity_changed=False,
                                          call_orbit_dissipation=False)

        # Now that those changes have been passed and the tides have been updated, we can make a single call to
        #    the orbit's dissipation_changed method
        self.dissipation_changed(world_instance)

    def dissipation_changed(self, world_signature: WorldSignatureType):
        """ Tidal dissipation has changed on the provided world. Make any necessary changes. """

        world_instance = self.world_signature_to_index(world_signature)
        log.debug(f'Dissipation changed for {world_instance} in {self}.')

        # LEFTOFF - Add dual dissipation calculations. de/dt and da/dt

    def add_star(self, star_world: StarWorld, is_tidal_host: bool = False, run_update: bool = True):
        """ Add a star to the orbit. This star may or may not be the tidal host.

        Stars that are not tidal hosts are only used only for insolation calculations.

        Parameters
        ----------
        star_world : StarWorld
            TidalPy star instance to be added to the orbit.
        is_tidal_host : bool = False
            If `True`, then the star will be added as a tidal world as well.
        run_update : bool = True
            If `True`, the orbit's `update_orbit` method will be called after the world has been added.
        """

        super().add_star(star_world, is_tidal_host=is_tidal_host, run_update=False)

        # If there are any tidal objects loaded into the orbit, calculate their insolation heating.
        okay_to_calculate_insolation = False
        if is_tidal_host:
            okay_to_calculate_insolation = True
        else:
            if self.tidal_host is not None:
                okay_to_calculate_insolation = True
        if okay_to_calculate_insolation:
            for tidal_world in self.tidal_objects:
                if tidal_world is self.star:
                    continue
                if self.get_stellar_distance(tidal_world) is not None:
                    self.calculate_insolation(tidal_world)

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
        super().add_tidal_world(tidal_world, is_tidal_host=is_tidal_host, run_update=False)

        # If the world is tidal then add it to the tidally active list
        if isinstance(tidal_world, all_tidal_world_types):
            # Check if the world can even be tidally active or not.
            if tidal_world.tides_on:
                # Check if tides are turned on or not inside world.
                if tidal_world not in self.tidally_active_worlds:
                    self._tidally_active_worlds.append(tidal_world)

        # Calculate insolation heating
        if self.star is not None:
            if tidal_world is not self.star:
                if self.get_stellar_distance(tidal_world) is not None:
                    self.calculate_insolation(tidal_world)

        # Initialize the tides model
        if is_tidal_host:
            # If the tidal host is added, try to initialize all world_types.
            for world in self.tidal_objects:
                if world is self.tidal_host:
                    # Need to have the host tide raiser loaded before we can initialize the tidal host's tides
                    if self.host_tide_raiser is not None and world.tides is not None:
                        world.tides.post_orbit_initialize()
                else:
                    if world.tides is not None:
                        world.tides.post_orbit_initialize()
        else:
            # For non-tidal host world_types, simply try to initialize the tides on only the provided world.
            if tidal_world.tides is not None and self.tidal_host is not None:
                tidal_world.tides.post_orbit_initialize()

        # Update orbit if needed
        if run_update:
            orbital_freq_changed = self.get_orbital_frequency(tidal_world) is not None
            eccentricity_changed = self.get_eccentricity(tidal_world) is not None
            self.orbit_changed(tidal_world, orbital_freq_changed=orbital_freq_changed,
                               eccentricity_changed=eccentricity_changed)

    def set_stellar_distance(self, world_signature: WorldSignatureType, distance: FloatArray):
        """ Set the orbital distance between a world of interest and the star (used for insolation calculations)

        If the tidal host is a star then this will simply wrap the world's semi-major axis setter. For a non-star host,
            then we assume that the world will share its stellar distance with its tidal host. For example, Io's
            solar flux is largely determined by Jupiter's orbit, not Io's orbit around Jupiter.

        Parameters
        ----------
        world_signature : WorldSignatureType
            A signature used to distinguish one tidal world from another. This could be its name,
                orbital location index, or the instance of an initialized TidalPy world.
        distance : FloatArray
            New stellar distance from the desired world to the host star [m].

        """

        super().set_stellar_distance(world_signature, distance)

        if self.star_host:
            # We do not need to recalculate insolation for all the world_types, only the one provided.
            self.calculate_insolation(world_signature)
        else:
            # If the stellar distance has changed and it is not a star host then we need to recalculate insolation
            #    for all world_types.
            for world in self.tidal_objects:
                self.calculate_insolation(world)

    def set_stellar_eccentricity(self, world_signature: WorldSignatureType, eccentricity: FloatArray):
        """ Set the orbital eccentricity between a world of interest and the star (used for insolation calculations)

        If the tidal host is a star then this will simply wrap the world's eccentricity setter. For a non-star host,
            then we assume that the world will share its stellar distance with its tidal host. For example, Io's
            solar flux is largely determined by Jupiter's orbit, not Io's orbit around Jupiter.

        Parameters
        ----------
        world_signature : WorldSignatureType
            A signature used to distinguish one tidal world from another. This could be its name,
                orbital location index, or the instance of an initialized TidalPy world.
        eccentricity : FloatArray
            New stellar eccentricity between the desired world to the host star [m].

        """

        super().set_stellar_distance(world_signature, eccentricity)

        if self.star_host:
            # We do not need to recalculate insolation for all the world_types, only the one provided.
            self.calculate_insolation(world_signature)
        else:
            # If the stellar distance has changed and it is not a star host then we need to recalculate insolation
            #    for all world_types.
            for world in self.tidal_objects:
                self.calculate_insolation(world)

    def calculate_insolation(self, world_signature: WorldSignatureType, set_insolation_in_world: bool = True):
        """ Calculate the insolation heating received by a world with the provided signature.

        The star-world separation (semi-major axis) and eccentricity are used to estimate the orbit-averaged
            insolation heating received at the surface of the target world. If the star is the host of the system then
            the target world's semi-major axis and eccentricity will be used.
            Otherwise the tidal host's parameters will be used.

        The actual method used to make the calculation is stored in the target world's equilibrium_insolation_func.
            It is set by the world's user-provided configuration---the only difference between the models is how
            they handle an eccentric orbit.

        Parameters
        ----------
        world_signature : WorldSignatureType
            A signature used to distinguish one tidal world from another. This could be its name,
                orbital location index, or the instance of an initialized TidalPy world.
        set_insolation_in_world : bool = True
            If `True`, method will make a call to the world's set_insolation method.

        """

        if self.star is None:
            insolation_heating = None
        else:
            # Get the world's instance
            orbit_index = self.world_signature_to_index(world_signature)
            world = self.tidal_objects[orbit_index]

            # Get its distance to the star (this will be the world's current semi-major axis if the tidal host is the star)
            #   Otherwise it will be the stellar distance of the tidal host to the star. This is all handled in the
            #   get_stellar_distance method.
            stellar_distance = self.get_stellar_distance(world_signature)
            stellar_eccentricity = self.get_stellar_eccentricity(world_signature)
            if stellar_eccentricity is None:
                # Stellar eccentricity not set. Assume it is zero
                log.warning(f'Calculating insolation heating for {world}, but stellar eccentricity is not set. '
                            f'Assuming zero.')
                stellar_eccentricity = np.zeros_like(stellar_distance)

            insolation_heating = \
                world.equilibrium_insolation_func(self.star.luminosity, stellar_distance, world.albedo, world.radius,
                                                  stellar_eccentricity)

            if set_insolation_in_world:
                world.set_insolation_heating(insolation_heating)

        return insolation_heating


    # # Initialized properties
    @property
    def tidally_active_worlds(self) -> List[TidalWorldType]:
        """ List of world_types stored in orbit that are tidally active. """
        return self._tidally_active_worlds

    @tidally_active_worlds.setter
    def tidally_active_worlds(self, value):
        raise InitiatedPropertyChangeError
