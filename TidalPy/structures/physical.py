import numpy as np

from TidalPy import debug_mode
from TidalPy.types import float_like, float_eps
from TidalPy.exceptions import (BadAttributeValueError, IncorrectAttributeType, ImproperAttributeHandling,
                                UnusualRealValueError, ParameterMissingError, TidalPyException)
from TidalPy.utilities.classes import ConfigHolder

from ..configurations import raise_on_changed_config
from .. import log


class ImproperAttributeChanged(TidalPyException):

    default_message = 'A pre-computed planet had a critical configuration change that requires a new instance to be created'

    def __init__(self, *args, force_raise: bool = False, **kwargs):

        # If no input is provided then the base exception will look at the class attribute 'default_message'
        #   and send that to sys.stderr

        if not force_raise and not raise_on_changed_config:
            log.warn(self.default_message)
        else:
            if args or kwargs:
                super().__init__(*args, **kwargs)
            else:
                super().__init__(self.default_message)


class PhysicalObjSpherical(ConfigHolder):

    """ PhysicalObj Class: contains attributes and functionality used for a few physical objects such as
    planets or layers

    Assumes spherical geometry
    """

    def __init__(self, config: dict, automate: bool = False):

        super().__init__(replacement_config=config, automate=automate)

        self.geometry_init = False

        # Required parameters that are set in the 'set_geometry' method
        self._mass = None
        self._radius = None
        self._thickness = None

        # parameters that are calculated
        self._radius_inner = None
        self._volume = None
        self._surface_area_outer = None
        self._surface_area_inner = None
        self._density_bulk = None
        self._moi = None
        self._moi_ideal = None
        self._moi_factor = None

        # State Variables
        self._time = None
        self._spin_freq = None
        self._orbital_freq = None

    def set_geometry(self, radius: float, mass: float, thickness: float = None):
        """ Sets and calculates object's physical parameters based on user provided input.

        :param radius:      <float> outer radius of object
        :param mass:        <float> mass of object
        :param thickness:   <float> thickness of object
        """

        if thickness is None:
            raise ParameterMissingError('Base class of PhysicalObjSpherical requires thickness to set geometry.')

        if debug_mode:
            for arg in [radius, thickness, mass]:
                if type(arg) not in float_like:
                    raise IncorrectAttributeType
                elif arg < 0:
                    raise BadAttributeValueError

        self._radius = radius
        self._thickness = thickness
        self._mass = mass

        # Update Physical Properties
        self._radius_inner = self.radius - self.thickness
        self._volume = (4. / 3.) * np.pi * (self.radius**3 - self.radius_inner**3)
        self._surface_area_outer = 4. * np.pi * self.radius**2
        self._surface_area_inner = 4. * np.pi * self.radius_inner**2
        self._density_bulk = self.mass / self.volume
        self._moi_ideal = (2. / 5.) * self.mass * (self.radius**5 - self.radius_inner**5) / \
                          (self.radius**3 - self.radius_inner**3)

        if self.moi is not None:
            self._moi_factor = self.moi / self.moi_ideal

        # Perform sanity checks
        if debug_mode:
            # Basic Checks
            if self.thickness > self.radius:
                raise BadAttributeValueError
            if self.radius_inner > self.radius:
                # Inner radius should be <= self.radius
                raise BadAttributeValueError
            elif abs(self.radius_inner - self.radius) < float_eps and self.thickness > float_eps:
                # Inner radius should be < self.radius is self.thickness != 0
                raise BadAttributeValueError

            # Realistic Value Checks
            if self.radius < 1. or self.radius > 1.e12:
                raise UnusualRealValueError
            if self.mass < 50. or self.mass > 1.e31:
                raise UnusualRealValueError
            if self.density_bulk < .3 or self.density_bulk > 1.e5:
                raise UnusualRealValueError

        self.geometry_init = True

    # State Variables
    @property
    def time(self):
        return self._time

    @time.setter
    def time(self, value: np.ndarray):
        if debug_mode:
            assert type(value) == np.ndarray
        self._time = value

    @property
    def spin_freq(self):
        return self._spin_freq

    @spin_freq.setter
    def spin_freq(self, value: np.ndarray):
        if debug_mode:
            assert type(value) == np.ndarray
        self._spin_freq = value

    @property
    def orbital_freq(self):
        return self._orbital_freq

    @orbital_freq.setter
    def orbital_freq(self, value: np.ndarray):
        if debug_mode:
            assert type(value) == np.ndarray
        self._orbital_freq = value

    # Primary Parameters
    @property
    def mass(self):
        return self._mass

    @mass.setter
    def mass(self, value):
        raise ImproperAttributeHandling('Mass must be set by the set_geometry method')

    @property
    def radius(self):
        return self._radius

    @radius.setter
    def radius(self, value):
        raise ImproperAttributeHandling('Radius must be set by the set_geometry method')

    @property
    def thickness(self):
        return self._thickness

    @thickness.setter
    def thickness(self, value):
        raise ImproperAttributeHandling('Thickness must be set by the set_geometry method')

    # Secondary Parameters
    @property
    def volume(self):
        return self._volume

    @volume.setter
    def volume(self, value):
        raise ImproperAttributeHandling('Volume is set by the set_geometry method')

    @property
    def radius_inner(self):
        return self._radius_inner

    @radius_inner.setter
    def radius_inner(self, value):
        raise ImproperAttributeHandling('Inner Radius is set by the set_geometry method')

    @property
    def surface_area_outer(self):
        return self._surface_area_outer

    @surface_area_outer.setter
    def surface_area_outer(self, value):
        raise ImproperAttributeHandling('Outer Surface Area is set by the set_geometry method')

    @property
    def surface_area_inner(self):
        return self._surface_area_inner

    @surface_area_inner.setter
    def surface_area_inner(self, value):
        raise ImproperAttributeHandling('Inner Surface Area is set by the set_geometry method')

    @property
    def density_bulk(self):
        return self._density_bulk

    @density_bulk.setter
    def density_bulk(self, value):
        raise ImproperAttributeHandling('Bulk Density is set by the set_geometry method')

    @property
    def moi_ideal(self):
        """ Ideal Moment of Inertia (spherical shell)"""
        return self._moi_ideal

    @moi_ideal.setter
    def moi_ideal(self, value):
        raise ImproperAttributeHandling('Ideal Moment of Inertia is set by the set_geometry method')

    @property
    def moi_factor(self):
        """ Moment of Inertia Factor (real moi / ideal moi) """
        return self._moi_factor

    @moi_factor.setter
    def moi_factor(self, value):
        raise ImproperAttributeHandling('Moment of Inertia Factor is set by the set_geometry method or '
                                        'by setting self.moi (the real moment of inertia)')

    @property
    def moi(self):
        """ Real Moment of Inertia (usually found via BurnMan for tidal planets)"""
        return self._moi

    @moi.setter
    def moi(self, value):

        # Moment of inertia could, in general, be calculated after an object is created so it may be set after both
        #   __init__ and set_geometry are called.
        if debug_mode:
            if type(value) not in float_like:
                raise IncorrectAttributeType
            elif value < 0.:
                raise BadAttributeValueError

        self._moi = value
        if self.moi_ideal is not None:
            self._moi_factor = self.moi / self.moi_ideal

    # Alias Attributes
    @property
    def surface_area(self):
        # Aliased with self.surface_area_outer
        return self.surface_area_outer

    @surface_area.setter
    def surface_area(self, value):
        self.surface_area_outer = value

    @property
    def M(self):
        # Aliased with self.mass
        return self.mass

    @M.setter
    def M(self, value):
        self.mass = value

    @property
    def V(self):
        # Aliased with self.volume
        return self.volume

    @V.setter
    def V(self, value):
        self.volume = value

    @property
    def R(self):
        # Aliased with self.radius
        return self.radius

    @R.setter
    def R(self, value):
        self.radius = value

    @property
    def dx(self):
        # Aliased with self.thickness
        return self.thickness

    @dx.setter
    def dx(self, value):
        self.thickness = value