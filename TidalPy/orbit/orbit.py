


class OrbitBase:


    def __init__(self, star, target_bodies: list, host = None):
        """

        :param star:
        :param host:
        :param target_bodies:
        """

        self.star = star
        if host is None:
            host = star
        self.host = host
        self.target_bodies = target_bodies

        # State variables
        self._eccentricities = list()
        self._inclinations = list()
        self._orbital_freqs = list()
        self._semi_major_axis = list()

        # Are target bodies orbiting the star or the host
        self.star_host = False
        if self.star is self.host:
            self.star_host = True

        # Determine the equilibrium temperatures on the tidal bodies
        for target_body in self.target_bodies:
            if self.star_host:
                seperation = target_body.semi_major_axis
            else:
                seperation = self.host.semi_major_axis
            target_body.temperature_equilibrium = calc_equilib_temp(seperation)
            # Store dummy values for the state variables
            self._eccentricities.append(None)
            self._inclinations.append(None)
            self._orbital_freqs.append(None)
            self._semi_major_axis.append(None)


        # Add a reference to the orbit to the planet(s)
        self.host._orbit = self
        self.star._orbit = self
        for t_i, target_body in enumerate(self.target_bodies):
            target_body._orbit = self
            target_body.orbit_num = t_i

