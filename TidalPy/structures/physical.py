""" physical.py - Physical object base class

This module contains the base python class for physical objects (layers, planets, stars, etc.).

"""

from typing import Union

import numpy as np

from .. import debug_mode, log
from ..constants import G
from ..exceptions import (BadAttributeValueError, ImproperPropertyHandling, IncorrectAttributeType,
                          UnusualRealValueError, MissingArgumentError, ImproperGeometryPropertyHandling)
from ..utilities.classes import ConfigHolder
from ..utilities.types import float_eps, float_like, NoneType


FloatNone = Union[NoneType, float]

class PhysicalObjSpherical(ConfigHolder):
    """ PhysicalObjSpherical Class contains attributes and functionality used for objects such as planets or layers
    that are spherical shell.

    Assumptions
    -----------
    Assumes spherical geometry.

    """

    def __init__(self, config: dict):

        super().__init__(replacement_config=config, store_py_info=True)

        # Configuration properties
        self._num_slices = None

        # Properties that are set in the 'set_geometry' method
        self._mass = None
        self._radius = None
        self._thickness = None
        self._mid_slice_index = None
        self.geometry_init = False
        self.moi_is_ideal = False

        # Properties that are calculated using the above parameters
        self._radius_inner = None
        self._radius_middle = None
        self._volume = None
        self._surface_area_outer = None
        self._surface_area_inner = None
        self._surface_area_middle = None
        self._density_bulk = None
        self._moi = None
        self._moi_ideal = None
        self._moi_factor = None
        self._beta_inner = None
        self._beta_middle = None
        self._beta_outer = None

        # Other properties than an object may or may not have depending on how it was initialized
        self._temperature_outer = None
        self._temperature_inner = None
        self._density_outer = None
        self._density_middle = None
        self._density_inner = None
        self._gravity_outer = None
        self._gravity_inner = None
        self._gravity_middle = None
        self._pressure_outer = None
        self._pressure_inner = None
        self._pressure_middle = None

        # Slice properties
        self._radii = None
        self._volume_slices = None
        self._sa_slices = None
        self._depths = None
        self._mass_slices = None
        self._mass_below_slices = None
        self._density_slices = None
        self._gravity_slices = None
        self._pressure_slices = None

        # Properties set by other constructors
        self._mass_below = None
        self._pressure_above = None

    def reinit(self, initial_init: bool = False, set_by_burnman: bool = False):
        """ Reinitialize the physical object by pulling in any potentially new configurations

        Parameters
        ----------
        initial_init : bool = False
            Set to `True` for the first time an instance is created.
        set_by_burnman : bool = False
            Set to `True` if a Burnman layer/world constructor is calling reinit
        """

        if initial_init:
            log.debug(f'First initialization called for {self}.')
        else:
            log.debug(f'Reinit called for {self}.')
            # If this is a reinit, then the state of the world should be cleared.
            self.clear_state()

        if not set_by_burnman:
            # Setup Geometry
            #    Pull out densities and pressures and convert them into constant value slices
            self._num_slices = self.config.get('slices', None)
            #    If user provided real moment of inertia, pull that out and calculate moi factor
            if self.config.get('moi', None) is not None:
                self.moi = self.config.get('moi', None)

        # Other reinit steps are set by child class' reinit methods.

    def set_geometry(self, radius: float, mass: float, thickness: float = None,
                     mass_below: float = 0., update_state_geometry: bool = True, build_slices: bool = True):
        """ Calculates and sets the object's physical parameters based on user provided input.

        Assumptions
        -----------
        Spherical Geometry

        Parameters
        ----------
        radius : float
            Outer radius of object [m]
        mass : float
            Mass of object [kg]
        thickness : float = None
            Thickness of the object [m]
        mass_below : float = 0.
            Mass below this object (only applicable for shell-like structures)
            Used in gravity and pressure calculations
        update_state_geometry : bool = True
            Update the class' state geometry
        build_slices : bool = True
            If True, method will attempt to calculate gravities, densities, etc. for each slice.

        """

        log.debug(f'Set geometry called for {self}')
        if thickness is None:
            # TODO: Why is the below error here?
            log.error('Base class of PhysicalObjSpherical requires thickness to set geometry.')
            raise MissingArgumentError('Base class of PhysicalObjSpherical requires thickness to set geometry.')

        if debug_mode:
            for arg in [radius, thickness, mass]:
                if type(arg) not in float_like:
                    raise IncorrectAttributeType
                elif arg < 0:
                    raise BadAttributeValueError

        # Set state properties (may not be called for burnman layers)
        if update_state_geometry:
            # Set the geometry state properties
            self._radius = radius
            self._thickness = thickness

            # Set the mass properties
            self._mass = mass
            self._mass_below = mass_below

            # Update physical properties
            self._radius_middle = self.radius - (self.thickness / 2.)
            self._radius_inner = self.radius - self.thickness
            self._volume = (4. / 3.) * np.pi * (self.radius**3 - self.radius_inner**3)
            self._surface_area_outer = 4. * np.pi * self.radius**2
            self._surface_area_middle = 4. * np.pi * self.radius_middle**2
            self._surface_area_inner = 4. * np.pi * self.radius_inner**2

            # Base models assume constant density throughout object
            self._density_bulk = self.mass / self.volume
            self._density_outer = self.density_bulk
            self._density_middle = self.density_bulk
            self._density_inner = self.density_bulk

            # Acceleration due to gravity is based on this layer/objects mass and any mass below it (if applicable)
            self._gravity_outer = G * (self.mass + self.mass_below) / self.radius**2
            half_volume_lower = (4. / 3.) * np.pi * (self.radius_middle**3 - self.radius_inner**3)
            half_mass_lower = half_volume_lower * self.density_bulk
            self._gravity_middle = G * (half_mass_lower + self.mass_below) / self.radius_middle**2
            if thickness is None or thickness == self.radius:
                self._gravity_inner = 0.
            else:
                self._gravity_inner = G * self.mass_below / self.radius_inner**2

            if build_slices and self.num_slices is not None:
                # Construct slices of the structure from the base to the top
                #    These are defined by their radius so the first most slice will have a radius != the radius_inner
                dx_per_slice = self.thickness / self.num_slices
                starting_radius = self.radius_inner + dx_per_slice
                self._radii = np.linspace(starting_radius, self.radius, self.num_slices, endpoint=True)
                self._depths = self.radius - self.radii
                self._sa_slices = 4. * np.pi * self.radii**2
                self._volume_slices = np.zeros_like(self._radii)
                self._volume_slices[0] = (4. / 3.) * np.pi * (self.radii[0]**3 - self.radius_inner**3)
                self._volume_slices[1:] = (4. / 3.) * np.pi * (self.radii[1:]**3 - self.radii[:-1]**3)

                # Base models assume constant density throughout object
                self._density_slices = self.density_bulk * np.ones_like(self.radii)
                self._mass_slices = self.density_slices * self.volume_slices

                # Mass below each slice is equal to slice masses + and mass below this physical object
                self._mass_below_slices = np.asarray(
                        [self.mass_below + sum(self.mass_slices[:i+1]) for i in range(self.num_slices)]
                )

                # Calculate gravity at the top of each slice
                self._gravity_slices = G * self.mass_below_slices / self.radii**2

        # Other properties that should be set for all physical methods
        self._beta_outer = self.gravity_outer * self.radius * self.density_outer
        self._beta_middle = self.gravity_middle * self.radius_middle * self.density_middle
        self._beta_inner = self.gravity_inner * self.radius_inner * self.density_inner
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

            # Slice Checks
            if build_slices and self.num_slices > 10:
                # If the num_slices < 10 then small errors may make these checks always fail.
                np.testing.assert_approx_equal(sum(self.mass_slices), self.mass)
                np.testing.assert_approx_equal(self.mass_below_slices[-1], self.mass_below + self.mass)
                np.testing.assert_approx_equal(self.gravity_slices[-1], self.gravity_outer)

        # Let the class know that the geometry has been set
        self.geometry_init = True

    def set_static_pressure(self, pressure_above: float = None, build_slices: bool = True,
                            called_from_burnman: bool = False):
        """ Sets the static pressure for the physical structure.

        `Static` here indicates that this is not a dynamic pressure used in many calculations. The static pressure can
            be used in place of the dynamic pressure, but that is not always the case.

        Parameters
        ----------
        pressure_above : float = None
            Pressure above this structure. If this is a layer, then it is the pressure at the base of the overlying
                layer. If it is the upper-most layer or a world, then it may be the surface pressure.
        build_slices : bool = True
            If `True`, method will find the pressure at each slice of the physical object.
        called_from_burnman : bool = False
            Set to `True` if called from a burnman layer/world.
        """

        log.debug(f'Setting up static pressure for {self}.')

        if called_from_burnman:
            log.debug('Static pressure method called from burnman - skipping.')
            return True

        if pressure_above is None:
            pressure_above = self.pressure_above

        if pressure_above is None:
            log.error(f'Not enough information to build static pressure for {self}.')
            raise MissingArgumentError

        else:
            # Calculate pressures from top down
            grav_radius_outer = self.gravity_outer * self.radius_outer
            grav_radius_middle = self.gravity_middle * self.radius_middle
            grav_radius_inner = self.gravity_inner * self.radius_inner
            self._pressure_outer = pressure_above
            self._pressure_middle = self.pressure_outer + \
                                    self.density_middle * (grav_radius_outer - grav_radius_middle)
            self._pressure_inner = self.pressure_middle + \
                                   self.density_inner * (grav_radius_middle - grav_radius_inner)

            if build_slices:
                # Calculate gravity * radius difference for each slice
                gravrad_diff_top = self.gravity_slices[1:] * self.radii[1:] - self.gravity_slices[:-1] * self.radii[:-1]
                # The bottom-most slice uses the radius and gravity at the bottom of the bottom-most slice
                gravrad_diff_bot = [self.gravity_slices[0] * self.radii[0] - self.gravity_inner * self.radius_inner]
                gravrad_diff = np.concatenate((gravrad_diff_bot, gravrad_diff_top))

                # Pressure pieces describe each slice's contribution to the pressure
                pressure_pieces = self.density_slices * gravrad_diff

                # The real pressure equals the pressure piece of a slice + any pressure above it
                self._pressure_slices = \
                    np.asarray(
                        [pressure_above + sum(pressure_pieces[i:]) for i in range(self.num_slices)]
                    )

                if debug_mode:
                    np.testing.assert_approx_equal(self.pressure_inner, self.pressure_slices[0])


    # # Independent state properties
    @property
    def mid_slice_index(self) -> int:
        """ The Slice Index Closest to (Radius - Thickness / 2) """
        return self._mid_slice_index

    @mid_slice_index.setter
    def mid_slice_index(self, value):
        raise ImproperGeometryPropertyHandling

    @property
    def mass(self) -> float:
        """ Physical Object's Mass [kg] """

        return self._mass

    @mass.setter
    def mass(self, value):
        raise ImproperGeometryPropertyHandling

    @property
    def radius(self) -> float:
        """ Physical Object's Outer Radius [m] """

        return self._radius

    @radius.setter
    def radius(self, value):
        raise ImproperGeometryPropertyHandling

    @property
    def thickness(self) -> float:
        """ Physical Object's Thickness [m] """

        return self._thickness

    @thickness.setter
    def thickness(self, value):
        raise ImproperGeometryPropertyHandling


    # Dependent state properties
    #    Geometry properties
    @property
    def volume(self) -> float:
        """ Physical Object's Volume [m^3] """

        return self._volume

    @volume.setter
    def volume(self, value):
        raise ImproperGeometryPropertyHandling

    @property
    def radius_inner(self) -> float:
        """ Physical Object's Inner Radius [m] """
        return self._radius_inner

    @radius_inner.setter
    def radius_inner(self, value):
        raise ImproperGeometryPropertyHandling

    @property
    def radius_middle(self) -> float:
        """ Physical Object's Middle Radius [m]

        "middle" is defined by the object's thickness / 2
        """

        return self._radius_middle

    @radius_middle.setter
    def radius_middle(self, value):
        raise ImproperGeometryPropertyHandling

    @property
    def surface_area_outer(self) -> float:
        """ Physical Object's Outer Surface Area [m-2] """

        return self._surface_area_outer

    @surface_area_outer.setter
    def surface_area_outer(self, value):
        raise ImproperGeometryPropertyHandling

    @property
    def surface_area_middle(self) -> float:
        """ Physical Object's Middle Surface Area [m-2]

        "middle" is defined by the object's thickness / 2
        """

        return self._surface_area_middle

    @surface_area_middle.setter
    def surface_area_middle(self, value):
        raise ImproperGeometryPropertyHandling

    @property
    def surface_area_inner(self) -> float:
        """ Physical Object's Inner Surface Area [m-2] """

        return self._surface_area_inner

    @surface_area_inner.setter
    def surface_area_inner(self, value):
        raise ImproperGeometryPropertyHandling

    @property
    def density_bulk(self) -> float:
        """ Physical Object's Bulk Density [kg m-3] """

        return self._density_bulk

    @density_bulk.setter
    def density_bulk(self, value):
        raise ImproperGeometryPropertyHandling

    @property
    def moi_ideal(self) -> float:
        """ Physical Object's Moment of Inertia (for an Ideal Sphere or Spherical Shell) [kg m^2]

        .. math:: C_{ \text{Ideal} } = (2 / 5) * M * (R^5 - R_{ \text{Inner} }^5) / (R^3 - R_{ \text{Inner} }^3)

        Assumptions
        -----------
        Spherical Geometry

        See Also
        --------
        PhysicalObjSpherical.moi
        PhysicalObjSpherical.moi_factor
        """

        return self._moi_ideal

    @moi_ideal.setter
    def moi_ideal(self, value):
        raise ImproperGeometryPropertyHandling

    @property
    def moi_factor(self) -> float:
        """ Physical Object's Moment of Inertia Factor [unitless]

        .. math:: C_{ \text{f} } = C / C_{ \text{Ideal} }

        `C` may either be a measured moment of inertia, or one that is calculate using a more rigorous
            method than moi_ideal.

        See Also
        --------
        PhysicalObjSpherical.moi_ideal
        PhysicalObjSpherical.moi_factor
        """

        return self._moi_factor

    @moi_factor.setter
    def moi_factor(self, value):
        raise ImproperGeometryPropertyHandling('Moment of Inertia Factor is set by the set_geometry method or '
                                               'by setting self.moi (the real moment of inertia)')

    @property
    def beta_middle(self) -> float:
        """ Physical Object's Beta Parameter - Middle of Object

        Beta is related to the object's radius, bulk density, and gravitational acceleration.

        .. math:: \beta = R \rho g

        See Also
        --------
        PhysicalObjSpherical.beta
        PhysicalObjSpherical.beta_inner
        """

        return self._beta_middle

    @beta_middle.setter
    def beta_middle(self, value):
        raise ImproperGeometryPropertyHandling

    @property
    def beta_inner(self) -> float:
        """ Physical Object's Beta Parameter - Inner side of Object

        Beta is related to the object's radius, bulk density, and gravitational acceleration.

        .. math:: \beta = R \rho g

        See Also
        --------
        PhysicalObjSpherical.beta_middle
        PhysicalObjSpherical.beta
        """

        return self._beta_inner

    @beta_inner.setter
    def beta_inner(self, value):
        raise ImproperGeometryPropertyHandling

    @property
    def beta_outer(self) -> float:
        """ Physical Object's Beta Parameter - at Surface

        Beta is related to the object's radius, bulk density, and gravitational acceleration.

        .. math:: \beta = R \rho g

        See Also
        --------
        PhysicalObjSpherical.beta_middle
        PhysicalObjSpherical.beta_inner
        """

        return self._beta_outer

    @beta_outer.setter
    def beta_outer(self, value):
        raise ImproperGeometryPropertyHandling

    @property
    def density_outer(self) -> float:
        """ Density at the Outer Edge / Surface of the Object [kg m3] """
        return self._density_outer

    @density_outer.setter
    def density_outer(self, value):
        raise ImproperGeometryPropertyHandling

    @property
    def density_middle(self) -> float:
        """ Density at the Middle of the Object [kg m3] """
        return self._density_middle

    @density_middle.setter
    def density_middle(self, value):
        raise ImproperGeometryPropertyHandling

    @property
    def density_inner(self) -> float:
        """ Density at the Inner Edge of the Object [kg m3] """
        return self._density_inner

    @density_inner.setter
    def density_inner(self, value):
        raise ImproperGeometryPropertyHandling

    @property
    def gravity_outer(self) -> float:
        """ Acceleration due to Gravity at the Outer Edge / Surface of the Object [m s-2] """
        return self._gravity_outer

    @gravity_outer.setter
    def gravity_outer(self, value):
        raise ImproperGeometryPropertyHandling

    @property
    def gravity_middle(self) -> float:
        """ Acceleration due to Gravity at the Outer Edge / Surface of the Object [m s-2] """
        return self._gravity_middle

    @gravity_middle.setter
    def gravity_middle(self, value):
        raise ImproperGeometryPropertyHandling

    @property
    def gravity_inner(self) -> float:
        """ Acceleration due to Gravity at the Inner Edge of the Object [m s-2] """
        return self._gravity_inner

    @gravity_inner.setter
    def gravity_inner(self, value):
        raise ImproperGeometryPropertyHandling

    @property
    def pressure_inner(self) -> float:
        """ Pressure at the Inner Edge of the Object [Pa] """
        return self._pressure_inner

    @pressure_inner.setter
    def pressure_inner(self, value):
        raise ImproperGeometryPropertyHandling('Pressures are calculated using the set_static_pressure method.')

    @property
    def pressure_middle(self) -> float:
        """ Pressure at the Middle of the Object [Pa] """
        return self._pressure_middle

    @pressure_middle.setter
    def pressure_middle(self, value):
        raise ImproperGeometryPropertyHandling('Pressures are calculated using the set_static_pressure method.')

    @property
    def pressure_outer(self) -> float:
        """ Pressure at the Outer Edge / Surface of the Object [Pa] """
        return self._pressure_outer

    @pressure_outer.setter
    def pressure_outer(self, value):
        raise ImproperGeometryPropertyHandling('Pressures are calculated using the set_static_pressure method.')

    @property
    def radii(self) -> np.ndarray:
        """ Physical Object's Radius Slices (bottom to top) [m] """
        return self._radii

    @radii.setter
    def radii(self, value):
        raise ImproperGeometryPropertyHandling

    @property
    def depths(self) -> np.ndarray:
        """ Physical Object's Depth Slices (bottom to top) [m] """
        return self._depths

    @depths.setter
    def depths(self, value):
        raise ImproperGeometryPropertyHandling

    @property
    def volume_slices(self) -> np.ndarray:
        """ Physical Object's Volume Slices (bottom to top) [m3] """
        return self._volume_slices

    @volume_slices.setter
    def volume_slices(self, value):
        raise ImproperGeometryPropertyHandling

    @property
    def sa_slices(self) -> np.ndarray:
        """ Physical Object's Surface Area Slices (bottom to top) [m2] """
        return self._sa_slices

    @sa_slices.setter
    def sa_slices(self, value):
        raise ImproperGeometryPropertyHandling

    @property
    def mass_slices(self) -> np.ndarray:
        """ Mass of each Slice within the Physical Object [kg] """
        return self._mass_slices

    @mass_slices.setter
    def mass_slices(self, value):
        raise ImproperGeometryPropertyHandling

    @property
    def mass_below_slices(self) -> np.ndarray:
        """ Mass Below each Slice of the Physical Object [kg] """
        return self._mass_below_slices

    @mass_below_slices.setter
    def mass_below_slices(self, value):
        raise ImproperGeometryPropertyHandling

    @property
    def pressure_slices(self) -> np.ndarray:
        """ Physical Object's Pressure Slices (bottom to top) [Pa] """
        return self._pressure_slices

    @pressure_slices.setter
    def pressure_slices(self, value):
        raise ImproperGeometryPropertyHandling

    @property
    def density_slices(self) -> np.ndarray:
        """ Physical Object's Density Slices (bottom to top) [kg m3] """
        return self._density_slices

    @density_slices.setter
    def density_slices(self, value):
        raise ImproperGeometryPropertyHandling

    @property
    def gravity_slices(self) -> np.ndarray:
        """ Physical Object's Acceleration due to Gravity Slices (bottom to top) [m s-2] """
        return self._gravity_slices

    @gravity_slices.setter
    def gravity_slices(self, value):
        raise ImproperGeometryPropertyHandling

    # # State properties set by user or other methods
    @property
    def moi(self) -> FloatNone:
        """ Physical Object's Moment of Inertia [kg m^2]

        This may either be a measured moment of inertia, or one that is calculate using a more rigorous
            method than moi_ideal. In TidalPy this is usually set by BurnMan calculations.

        See Also
        --------
        PhysicalObjSpherical.moi_ideal
        PhysicalObjSpherical.moi_factor
        """

        if self._moi is None:
            if self.moi_ideal is not None:
                self.moi_is_ideal = False
                return self.moi_ideal
            else:
                return None
        else:
            self.moi_is_ideal = True
            return self._moi

    @moi.setter
    def moi(self, value):

        # Moment of inertia could, in general, be calculated after an object is created so it may be set after both
        #   __init__ and set_geometry are called.
        if debug_mode:
            if value is None:
                pass
            elif type(value) not in float_like:
                raise IncorrectAttributeType
            elif value < 0.:
                raise BadAttributeValueError

        self._moi = value
        if self.moi_ideal is not None:
            self._moi_factor = self.moi / self.moi_ideal

    @property
    def temperature_outer(self) -> float:
        """ Temperature at the Outer Edge / Surface of the Object [K] """
        return self._temperature_outer

    @temperature_outer.setter
    def temperature_outer(self, value):
        self._temperature_outer = value

    @property
    def temperature_inner(self) -> float:
        """ Temperature at the Inner Edge of the Object [K] """
        return self._temperature_inner

    @temperature_inner.setter
    def temperature_inner(self, value):
        self._temperature_inner = value

    @property
    def pressure_above(self) -> float:
        """ Pressure above this structure (surface pressure, or pressure from overlying layers) [Pa] """
        return self._pressure_above

    @pressure_above.setter
    def pressure_above(self, value):
        if debug_mode:
            assert type(value) == float

        self._pressure_above = value
        self.set_static_pressure()

    @property
    def mass_below(self) -> float:
        """ Mass below this structure (for worlds: this is zero; for layers: mass of the layer below) [kg] """
        return self._mass_below

    @mass_below.setter
    def mass_below(self, value):
        if debug_mode:
            assert type(value) == float

        self._mass_below = value
        self.set_static_pressure()

    # # Configuration properties
    @property
    def num_slices(self) -> int:
        """ Number of Slices cut into the Physical Object

        These are used for the calculation of pressures, gravities, etc. so that more accurate values can be found for
            e.g., pressure_inner, gravity_inner, etc.
        """

        return self._num_slices

    @num_slices.setter
    def num_slices(self, value):
        raise ImproperPropertyHandling('Number of slices must be set at initial construction.')

    # # Aliased properties
    # Global properties
    @property
    def gravity_surf(self) -> float:
        """ Alias of PhysicalObjSpherical.gravity_outer """
        return self.gravity_outer

    @gravity_surf.setter
    def gravity_surf(self, value):
        self.gravity_outer = value

    @property
    def beta(self) -> float:
        """ Alias of PhysicalObjSpherical.beta_outer

        This is the beta used in tidal calculations

        Assumptions
        -----------
        We have assumed that the beta generally needed is the beta at the surface of a layer/world.
        """

        # TODO: Should this be surface (as is used in Efroimsky work and others) or beta_middle?
        return self.beta_outer

    @beta.setter
    def beta(self, value):
        self.beta_outer = value

    @property
    def gravity_surface(self):
        """ Alias of PhysicalObjSpherical.gravity_outer [m s-2] """
        return self.gravity_outer

    @gravity_surface.setter
    def gravity_surface(self, value):
        self.gravity_outer = value

    @property
    def surface_area(self):
        """ Alias of PhysicalObjSpherical.surface_area_outer [m2] """
        return self.surface_area_outer

    @surface_area.setter
    def surface_area(self, value):
        self.surface_area_outer = value

    @property
    def radius_outer(self):
        """ Alias of PhysicalObjSpherical.radius [m] """
        return self.radius

    @radius_outer.setter
    def radius_outer(self, value):
        self.radius = value

    # Slice properties
    @property
    def volumes(self):
        """ Alias of PhysicalObjSpherical.volume_slices [m3] """
        return self.volume_slices

    @volumes.setter
    def volumes(self, value):
        self.volume_slices = value

    @property
    def masses(self):
        """ Alias of PhysicalObjSpherical.mass_slices [kg] """
        return self.mass_slices

    @masses.setter
    def masses(self, value):
        self.mass_slices = value

    @property
    def surface_areas(self):
        """ Alias of PhysicalObjSpherical.sa_slices [m2] """
        return self.sa_slices

    @sa_slices.setter
    def sa_slices(self, value):
        self.sa_slices = value

    @property
    def densities(self):
        """ Alias of PhysicalObjSpherical.density_slices [kg m-3] """
        return self.density_slices

    @densities.setter
    def densities(self, value):
        self.density_slices = value

    @property
    def gravities(self):
        """ Alias of PhysicalObjSpherical.gravity_slices [m s-2]"""
        return self.gravity_slices

    @gravities.setter
    def gravities(self, value):
        self.gravity_slices = value

    @property
    def pressures(self):
        """ Alias of PhysicalObjSpherical.pressure_slices [Pa] """
        return self.pressure_slices

    @pressures.setter
    def pressures(self, value):
        self.pressure_slices = value

    @property
    def M(self):
        """ Alias of PhysicalObjSpherical.mass [kg] """
        return self.mass

    @M.setter
    def M(self, value):
        self.mass = value

    @property
    def V(self):
        """ Alias of PhysicalObjSpherical.volume [m3] """
        return self.volume

    @V.setter
    def V(self, value):
        self.volume = value

    @property
    def R(self):
        """ Alias of PhysicalObjSpherical.radius [m] """
        return self.radius

    @R.setter
    def R(self, value):
        self.radius = value

    @property
    def dx(self):
        """ Alias of PhysicalObjSpherical.thickness [m] """
        return self.thickness

    @dx.setter
    def dx(self, value):
        self.thickness = value

    @property
    def g(self):
        """ Alias of PhysicalObjSpherical.gravity_outer [m s-2] """
        return self.gravity_outer

    @g.setter
    def g(self, value):
        self.gravity_outer = value
