import numpy as np

from .. import debug_mode
from ..constants import G
from ..exceptions import (BadAttributeValueError, ImproperPropertyHandling, IncorrectAttributeType,
                          UnusualRealValueError, ParameterMissingError)
from TidalPy.utilities.types import float_eps, float_like
from ..utilities.classes import ConfigHolder


class PhysicalObjSpherical(ConfigHolder):
    """ PhysicalObj Class contains attributes and functionality used for objects such as planets or
        spherical-shell layers

    Assumes spherical geometry
    """

    def __init__(self, config: dict):

        super().__init__(replacement_config=config)

        self.pyname = f'{self.__class__}'
        self.geometry_init = False

        # Properties that are set in the 'set_geometry' method
        self._mass = None
        self._radius = None
        self._thickness = None

        # Properties that are calculated using the above parameters
        self._radius_inner = None
        self._volume = None
        self._surface_area_outer = None
        self._surface_area_inner = None
        self._density_bulk = None
        self._moi = None
        self._moi_ideal = None
        self._moi_factor = None
        self._gravity_surf = None
        self._beta = None

    def reinit(self):
        """ Reinitialize the physical layer by pulling in any potentially new config changes.
        """

        # Setup by child classes
        pass

    def config_update(self):
        """ Any changes to the configuration of an object should be followed by a call to its config_update"""

        pass

    def set_geometry(self, radius: float, mass: float, thickness: float = None):
        """ Sets and calculates object's physical parameters based on user provided input.

        Parameters
        ----------
        radius : float
            Outer radius of object [m]
        mass : float
            Mass of object [kg]
        thickness : float = None
            Thickness of the object [m]
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
        self._gravity_surf = G * self.mass / self.radius**2
        self._radius_inner = self.radius - self.thickness
        self._volume = (4. / 3.) * np.pi * (self.radius**3 - self.radius_inner**3)
        self._surface_area_outer = 4. * np.pi * self.radius**2
        self._surface_area_inner = 4. * np.pi * self.radius_inner**2
        self._density_bulk = self.mass / self.volume
        self._moi_ideal = (2. / 5.) * self.mass * (self.radius**5 - self.radius_inner**5) / \
                          (self.radius**3 - self.radius_inner**3)
        self._beta = self.gravity_surf * self.radius * self.density_bulk

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
            elif abs(self.radius_inner - self.radius) < float_eps < self.thickness:
                # Inner radius should be < self.radius if self.thickness != 0
                raise BadAttributeValueError

            # Realistic Value Checks
            if self.radius < 1. or self.radius > 1.e12:
                raise UnusualRealValueError
            if self.mass < 50. or self.mass > 1.e31:
                raise UnusualRealValueError
            if self.density_bulk < .3 or self.density_bulk > 1.e5:
                raise UnusualRealValueError

        self.geometry_init = True

    # Independent state properties
    @property
    def mass(self) -> float:
        return self._mass

    @mass.setter
    def mass(self, value):
        raise ImproperPropertyHandling('Mass must be set by the set_geometry method')

    @property
    def radius(self) -> float:
        return self._radius

    @radius.setter
    def radius(self, value):
        raise ImproperPropertyHandling('Radius must be set by the set_geometry method')

    @property
    def thickness(self) -> float:
        return self._thickness

    @thickness.setter
    def thickness(self, value):
        raise ImproperPropertyHandling('Thickness must be set by the set_geometry method')

    # Dependent state properties
    @property
    def volume(self) -> float:
        return self._volume

    @volume.setter
    def volume(self, value):
        raise ImproperPropertyHandling('Volume is set by the set_geometry method')

    @property
    def radius_inner(self) -> float:
        return self._radius_inner

    @radius_inner.setter
    def radius_inner(self, value):
        raise ImproperPropertyHandling('Inner Radius is set by the set_geometry method')

    @property
    def surface_area_outer(self) -> float:
        return self._surface_area_outer

    @surface_area_outer.setter
    def surface_area_outer(self, value):
        raise ImproperPropertyHandling('Outer Surface Area is set by the set_geometry method')

    @property
    def surface_area_inner(self) -> float:
        return self._surface_area_inner

    @surface_area_inner.setter
    def surface_area_inner(self, value):
        raise ImproperPropertyHandling('Inner Surface Area is set by the set_geometry method')

    @property
    def density_bulk(self) -> float:
        return self._density_bulk

    @density_bulk.setter
    def density_bulk(self, value):
        raise ImproperPropertyHandling('Bulk Density is set by the set_geometry method')

    @property
    def gravity_surf(self) -> float:
        """ Surface Gravity of a sphere"""
        return self._gravity_surf

    @gravity_surf.setter
    def gravity_surf(self, value):
        raise ImproperPropertyHandling('Surface Gravity is set by the set_geometry method')

    @property
    def moi_ideal(self) -> float:
        """ Ideal Moment of Inertia (spherical shell)"""
        return self._moi_ideal

    @moi_ideal.setter
    def moi_ideal(self, value):
        raise ImproperPropertyHandling('Ideal Moment of Inertia is set by the set_geometry method')

    @property
    def moi_factor(self) -> float:
        """ Moment of Inertia Factor (real moi / ideal moi) """
        return self._moi_factor

    @moi_factor.setter
    def moi_factor(self, value):
        raise ImproperPropertyHandling('Moment of Inertia Factor is set by the set_geometry method or '
                                        'by setting self.moi (the real moment of inertia)')

    @property
    def beta(self) -> float:
        """ Beta Parameter (R * rho * g) """
        return self._beta

    @beta.setter
    def beta(self, value):
        raise ImproperPropertyHandling('Beta Parameter is set by the set_geometry method')

    @property
    def moi(self) -> float:
        """ Real Moment of Inertia (usually found via BurnMan for tidal planets)"""
        return self._moi

    # Properties set by user or other classes
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

    # Alias properties
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
