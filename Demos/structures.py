from typing import List, Tuple

import numpy as np

from TidalPy.utilities.conversions import days2rads, orbital_motion2semi_a, semi_a2orbital_motion
from TidalPy.rheology.base import RheologyModelBase
from TidalPy.rheology import Maxwell, Andrade, Sundberg, Elastic, Newton

# Build layer data storage structures
class LayerData:
    name: str
    radius: float
    layer_type: str
    is_static: bool
    is_incompressible: bool
    density: float
    shear_modulus: float
    bulk_modulus: float
    shear_viscosity: float
    bulk_viscosity: float
    shear_rheology_model: RheologyModelBase
    bulk_rheology_model: RheologyModelBase

    def __init__(
            self,
            radius: float,
            name: str = None,
            density: float = None,
            shear_modulus: float = None,
            bulk_modulus: float = None,
            shear_viscosity: float = None,
            bulk_viscosity: float = None):
        self.radius = radius
        if density is not None:
            self.density = density
        if shear_modulus is not None:
            self.shear_modulus = shear_modulus
        if bulk_modulus is not None:
            self.bulk_modulus = bulk_modulus
        if shear_viscosity is not None:
            self.shear_viscosity = shear_viscosity
        if bulk_viscosity is not None:
            self.bulk_viscosity = bulk_viscosity
        if name is not None:
            self.name = name

# Set default values for code
class CoreLayer(LayerData):
    name = "Core"
    layer_type = "solid"
    is_static = False
    is_incompressible = False
    shear_rheology_model = Maxwell()
    bulk_rheology_model = Elastic()

    # TODO: Find better approximations for the below.
    density = 20000.0
    shear_modulus = 500.0e9
    bulk_modulus = 1500.0e9
    shear_viscosity = 1.0e30
    bulk_viscosity = 1.0e30

class MantleLayer(LayerData):
    name = "Mantle"
    layer_type = "solid"
    is_static = False
    is_incompressible = False
    shear_rheology_model = Maxwell()
    bulk_rheology_model = Elastic()

    # TODO: Find better approximations for the below.
    density = 10000.0
    shear_modulus = 350.0e9
    bulk_modulus = 800.0e9
    shear_viscosity = 1.0e28
    bulk_viscosity = 1.0e30

class GasLayer(LayerData):
    name = "Gas Envelope"
    layer_type = "liquid"
    is_static = True
    is_incompressible = False
    shear_rheology_model = Newton()
    bulk_rheology_model = Elastic()

    # TODO: Find better approximations for the below.
    density = 800.0
    shear_modulus = 0.0
    bulk_modulus = 500.0e9
    shear_viscosity = 1.0e3
    bulk_viscosity = 1.0e30

# Planet data storage structures
class Planet:
    name: str
    radius: float
    _orbital_period: float  # in [days]
    _orbital_freq: float    # in [rad s-1]
    _semi_major_axis: float # in [m]
    spin_period: float      # in [days]
    host_mass: float        # in [kg]
    eccentricity: float
    obliquity: float

    layers: List[LayerData]

    def __init__(
            self,
            name: str, 
            radius: float,
            layers: List[LayerData],
            host_mass: float,
            eccentricity: float = 0.0,
            obliquity: float = 0.0,
            semi_major_axis: float = None,
            orbital_period: float = None,
            spin_period: float = None):

        self.name = name
        self.radius = radius
        self.layers = layers
        self.host_mass = host_mass
        self.eccentricity = eccentricity
        self.obliquity = obliquity

        if semi_major_axis is not None:
            self.semi_major_axis = semi_major_axis
        if orbital_period is not None:
            self.orbital_period = orbital_period

        if spin_period is None:
            # Assume tidally locked.
            self.spin_period = orbital_period
        else:
            self.spin_period = spin_period

    @property
    def spin_freq(self):
        if self.spin_period:
            return 2.0 * np.pi / (self.spin_period * 86400.0)
        else:
            return None

    @spin_freq.setter
    def spin_freq(self, value):
        raise AttributeError()
    
    @property
    def orbital_freq(self):
        return self._orbital_freq

    @orbital_freq.setter
    def orbital_freq(self, value):
        self._orbital_freq = value
        self._orbital_period = (2.0 * np.pi / self._orbital_freq) / 86400.0
        self._semi_major_axis = orbital_motion2semi_a(self._orbital_freq, self.host_mass)
    
    @property
    def orbital_period(self):
        return self._orbital_period
    
    @orbital_period.setter
    def orbital_period(self, value):
        self._orbital_period = value
        self._orbital_freq = 2.0 * np.pi / (self._orbital_period * 86400.0)
        self._semi_major_axis = orbital_motion2semi_a(self.orbital_freq, self.host_mass)

    @property
    def semi_major_axis(self):
        return self._semi_major_axis
    
    @semi_major_axis.setter
    def semi_major_axis(self, value):
        self._semi_major_axis = value
        self._orbital_freq = semi_a2orbital_motion(self._semi_major_axis, self.host_mass)
        self._orbital_period = 2.0 * np.pi / (self._orbital_freq * 86400.0)