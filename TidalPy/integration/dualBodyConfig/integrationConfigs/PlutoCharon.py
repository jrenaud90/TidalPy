import numpy as np
from scipy.constants import R

from ....radiogenics.defaults import standard_isotope_input
# Planet and Layer Radii are from Nimmo+2017

pluto = {
    'name':                'Pluto',
    'radius':              1188.3e3,
    'mass':                1.3895506e22,
    'constant_obliquity':  0.,
    'surface_temperature': 40.,
    'modern_semi_major_axis': 19.596e6,
    'tides_on': True,
    'moi':5.65e33,
    'layers':              {
        'Core':      {
            # Pre-calculated Layer Structure
            'radius_upper':                      853.2e3,
            'radius_lower':                      0.,
            'mass':                              3400. * (4. / 3.) * np.pi * 853.2e3**3,
            # Material Parameters
            'material_density':                  3400.,
            'thermal_expansion':                 5.0e-5,
            'thermal_conductivity':              4.0,
            'specific_heat':                     1225.5,
            'latent_heat':                       0.,
            'static_shear':                      50.e9,
            'solidus_temperature':               1800.,
            'liquidus_temperature':              2000.,
            # Convection Parameters
            'convection_alpha':                  1.0,
            'convection_beta':                   1. / 3.,
            'critical_rayleigh':                 1000.,
            # Solid Viscosity Model
            'viscosity_model':                   'reference',
            'viscosity_input':                   (1.5e14, 1800., 370.e3, 0.),
            # Partial Melt Model
            'partial_melt_model':                'off',
            'partial_melt_input':                tuple(),
            # Rheology / Complex Compliance Model
            'rheology_model':                    'andrade',
            'rheology_input':                    (.3, 1.),
            'tidal_scale':                       1.,
            'force_tides_off':                   False,
            # Other Model Inputs
            'growth_layer':                      False,
            'constant_ocean_temperature':        0.,
            'constant_viscoelastic_temperature': 0.,
            # Radiogenic Model
            'use_radiogenics':                   True,
            'radiogenic_model':                  'isotope',
            'radiogenic_input':                  standard_isotope_input

        },
        'Icy Shell': {
            # Pre-calculated Layer Structure
            'radius_upper':                      1189.9e3,
            'radius_lower':                      853.2e3,
            'mass':                              950. * (4. / 3.) * np.pi * (1189.9e3**3 - 853.2e3**3),
            # Material Parameters
            'material_density':                  950.,
            'thermal_expansion':                 1.56e-4,
            'thermal_conductivity':              2.27,
            'specific_heat':                     1925.,
            'latent_heat':                       284000.0,
            'static_shear':                      3.3e9,
            'solidus_temperature':               260.,
            'liquidus_temperature':              273.15,
            # Convection Parameters
            'convection_alpha':                  1.0,
            'convection_beta':                   1. / 3.,
            'critical_rayleigh':                 900.,
            # Solid Viscosity Model
            'viscosity_model':                   'reference',
            'viscosity_input':                   (5.e13, 260., 58.4e3, 0.),
            # Partial Melt Model
            'partial_melt_model':                'off',
            'partial_melt_input':                tuple(),
            # Rheology / Complex Compliance Model
            'rheology_model':                    'andrade',
            'rheology_input':                    (.3, 1.),
            'tidal_scale':                       1.,
            'force_tides_off':                   False,
            # Other Model Inputs
            'growth_layer':                      True,
            'constant_ocean_temperature':        273.15,
            'constant_viscoelastic_temperature': 273.,
            'viscoelastic_top_temperature':      273.15 / ((np.log(10.) / (58.4e3 / (R * 273.15))) + 1.),
            # Radiogenic Model
            'use_radiogenics':                   False,
            'radiogenic_model':                  'isotope',
            'radiogenic_input':                  standard_isotope_input
        }
    }
}

charon = {
    'name':                'Charon',
    'radius':              606000.0,
    'mass':                1.670299e21,
    'constant_obliquity':  0.,
    'surface_temperature': 40.,
    'modern_semi_major_axis': 19.596e6,
    'tides_on': True,
    'moi':1.77e32,
    'layers':              {
        'Core':      {
            # Pre-calculated Layer Structure
            'radius_upper':                      408.7e3,
            'radius_lower':                      0.,
            'mass':                              3400. * (4. / 3.) * np.pi * (408.7e3**3),
            # Material Parameters
            'material_density':                  3400.,
            'thermal_expansion':                 5.0e-5,
            'thermal_conductivity':              4.0,
            'specific_heat':                     1225.5,
            'latent_heat':                       0.,
            'static_shear':                      50.e9,
            'solidus_temperature':               1800.,
            'liquidus_temperature':              2000.,
            # Convection Parameters
            'convection_alpha':                  1.0,
            'convection_beta':                   1. / 3.,
            'critical_rayleigh':                 1000.,
            # Solid Viscosity Model
            'viscosity_model':                   'reference',
            'viscosity_input':                   (1.5e14, 1800., 370.e3, 0.),
            # Partial Melt Model
            'partial_melt_model':                'off',
            'partial_melt_input':                tuple(),
            # Rheology / Complex Compliance Model
            'rheology_model':                    'andrade',
            'rheology_input':                    (.3, 1.),
            'tidal_scale':                       1.,
            'force_tides_off':                   False,
            # Other Model Inputs
            'growth_layer':                      False,
            'constant_ocean_temperature':        0.,
            'constant_viscoelastic_temperature': 0.,
            # Radiogenic Model
            'use_radiogenics':                   True,
            'radiogenic_model':                  'isotope',
            'radiogenic_input':                  standard_isotope_input
        },
        'Icy Shell': {
            # Pre-calculated Layer Structure
            'radius_upper':                      606000.0,
            'radius_lower':                      408.7e3,
            'mass':                              950. * (4. / 3.) * np.pi * (606000.0**3 - 408.7e3**3),
            # Material Parameters
            'material_density':                  950.,
            'thermal_expansion':                 1.56e-4,
            'thermal_conductivity':              2.27,
            'specific_heat':                     1925.,
            'latent_heat':                       284000.0,
            'static_shear':                      3.3e9,
            'solidus_temperature':               260.,
            'liquidus_temperature':              273.15,
            # Convection Parameters
            'convection_alpha':                  1.0,
            'convection_beta':                   1. / 3.,
            'critical_rayleigh':                 900.,
            # Solid Viscosity Model
            'viscosity_model':                   'reference',
            'viscosity_input':                   (5.e13, 260., 58.4e3, 0.),
            # Partial Melt Model
            'partial_melt_model':                'off',
            'partial_melt_input':                tuple(),
            # Rheology / Complex Compliance Model
            'rheology_model':                    'andrade',
            'rheology_input':                    (.3, 1.),
            'tidal_scale':                       1.,
            'force_tides_off':                   False,
            # Other Model Inputs
            'growth_layer':                      True,
            'constant_ocean_temperature':        273.15,
            'constant_viscoelastic_temperature': 273.,
            'viscoelastic_top_temperature':      273.15 / ((np.log(10.) / (58.4e3 / (R * 273.15))) + 1.),
            # Radiogenic Model
            'use_radiogenics':                   False,
            'radiogenic_model':                  'isotope',
            'radiogenic_input':                  standard_isotope_input
        }
    }
}
