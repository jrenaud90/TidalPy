from TidalPy import debug_mode
import numpy as np
from typing import List, Union
from scipy.constants import G

from TidalPy.exceptions import ImproperAttributeHandling


class OrbitBase:


    def __init__(self, star, host = None, target_bodies: list = None):
        """

        :param star:
        :param host:
        :param target_bodies:
        """

        self.star = star
        if host is None:
            host = star
        self.host = host
        if target_bodies is None:
            self.target_bodies = list()
        else:
            self.target_bodies = target_bodies
        self.all_objects = [self.star, self.host] + self.target_bodies

        # State orbital variables (must be at least 2: for the star and the host)
        self._eccentricities = [None, None]  # type: List[Union[None, np.ndarray]]
        self._inclinations = [None, None]    # type: List[Union[None, np.ndarray]]
        self._orbital_freqs = [None, None]   # type: List[Union[None, np.ndarray]]
        self._semi_major_axis = [None, None] # type: List[Union[None, np.ndarray]]

        # Are target bodies orbiting the star or the host
        self.star_host = False
        if self.star is self.host:
            self.star_host = True

        # Determine various constants that depend upon other objects properties
        for target_body in self.target_bodies:

            # Equilibrium temperature
            # TODO
            # if self.star_host:
            #     separation = target_body.semi_major_axis
            # else:
            #     separation = self.host.semi_major_axis
            # target_body.temperature_equilibrium = calc_equilib_temp(separation)

            # Inflated Tidal Susceptibility
            target_body.tidal_susceptibility_inflated = (3. / 2.) * G * self.host.mass**2 * target_body.radius**5

            # Store dummy values for the state variables
            self._eccentricities.append(None)
            self._inclinations.append(None)
            self._orbital_freqs.append(None)
            self._semi_major_axis.append(None)

        # Add a reference to the orbit to the planet(s)
        self.host._orbit = self
        self.star._orbit = self
        self.host.orbit_location = 0
        self.host.orbit_location = 1
        for t_i, target_body in enumerate(self.target_bodies):
            target_body._orbit = self
            target_body.orbit_location = t_i + 2 # The +2 is due to the star==0 and host==1

    # Functionality to access and change the orbit's state variables

    def set_state(self, planet_loc: int, orbital_freq: np.ndarray = None, semi_major_axis: np.ndarray = None,
                  eccentricity: np.ndarray = None, inclination: np.ndarray = None, set_by_planet: bool = False):
        """ Set multiple orbital parameters at once, this reduces the number of calls to planet.orbit_change

        :param planet_loc:      <int> Orbital location of object's whose orbit is changing. Note: star=0, host=1
        :param orbital_freq:    <ndarray> (Optional, exclusive w/ semi_a)
        :param semi_major_axis: <ndarray> (Optional, exclusive w/ orb_freq)
        :param eccentricity:    <ndarray> (Optional)
        :param inclination:     <ndarray> (Optional)
        :param set_by_planet:   <bool> (=False) If true then it will not try to update other orbit-dependent variables
        """

        # The "set_by_planet" are set to True below even if the "set_by_planet" is False above because the function
        #    Will do a single call to orbit_update at the end.
        if orbital_freq is not None:
            if semi_major_axis is not None:
                raise ImproperAttributeHandling('Only set orbital frequency or semi-major axis, not both.')
            if type(orbital_freq) != np.ndarray:
                orbital_freq = np.asarray(orbital_freq)
            self.set_orbital_freq(planet_loc, orbital_freq, set_by_planet=True)
        if semi_major_axis is not None:
            if type(semi_major_axis) != np.ndarray:
                semi_major_axis = np.asarray(semi_major_axis)
            self.set_semi_major_axis(planet_loc, semi_major_axis, set_by_planet=True)
        if eccentricity is not None:
            if type(eccentricity) != np.ndarray:
                eccentricity = np.asarray(eccentricity)
            self.set_eccentricity(planet_loc, eccentricity, set_by_planet=True)
        if inclination is not None:
            if type(inclination) != np.ndarray:
                inclination = np.asarray(inclination)
            self.set_inclination(planet_loc, eccentricity, set_by_planet=True)

        if not set_by_planet:
            # Need to tell the planet to update any orbital-dependent state
            self.all_objects[planet_loc].orbit_update()


    def get_eccentricity(self, planet_loc: int):

        return self._eccentricities[planet_loc]

    def set_eccentricity(self, planet_loc: int, new_eccentricity: np.ndarray, set_by_planet: bool = False):

        if type(new_eccentricity) != np.ndarray:
            new_orbital_freq = np.asarray(new_eccentricity)

        self._eccentricities[planet_loc] = new_eccentricity

        if not set_by_planet:
            # Need to tell the planet to update any orbital-dependent state
            self.all_objects[planet_loc].orbit_update()

    def get_inclination(self, planet_loc: int):

        return self._inclinations[planet_loc]

    def set_inclination(self, planet_loc: int, new_inclination: np.ndarray, set_by_planet: bool = False):

        if type(new_inclination) != np.ndarray:
            new_orbital_freq = np.asarray(new_inclination)

        self._inclinations[planet_loc] = new_inclination

        if not set_by_planet:
            # Need to tell the planet to update any orbital-dependent state
            self.all_objects[planet_loc].orbit_update()

    def get_orbital_freq(self, planet_loc: int):

        return self._orbital_freqs[planet_loc]

    def set_orbital_freq(self, planet_loc: int, new_orbital_freq: np.ndarray, set_by_planet: bool = False):

        if type(new_orbital_freq) != np.ndarray:
            new_orbital_freq = np.asarray(new_orbital_freq)

        self._orbital_freqs[planet_loc] = new_orbital_freq

        # Changing the orbital frequency also changes the semi-major axis. Update via Kepler Laws
        if planet_loc not in [0, 1]:
            self._semi_major_axis[planet_loc] = np.cbrt(G * (self.all_objects[planet_loc].mass * self.host.mass) /
                                                        new_orbital_freq**2)
        if not set_by_planet:
            # Need to tell the planet to update any orbital-dependent state
            self.all_objects[planet_loc].orbit_update()

    def get_semi_major_axis(self, planet_loc: int):

        return self._semi_major_axis[planet_loc]

    def set_semi_major_axis(self, planet_loc: int, new_semi_major_axis: np.ndarray, set_by_planet: bool = False):

        if type(new_semi_major_axis) != np.ndarray:
            new_semi_major_axis = np.asarray(new_semi_major_axis)

        self._semi_major_axis[planet_loc] = new_semi_major_axis
        # Changing the orbital semi-major axis also changes the orbital frequency. Update via Kepler Laws
        if planet_loc not in [0, 1]:
            self._orbital_freqs[planet_loc] = np.sqrt(G * (self.all_objects[planet_loc].mass * self.host.mass) /
                                                        new_semi_major_axis**2)

        # Changing the orbital frequency also changes the semi-major axis. Update via Kepler Laws

        if not set_by_planet:
            # Need to tell the planet to update any orbital-dependent state
            self.all_objects[planet_loc].orbit_change()
