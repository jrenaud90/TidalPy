import numpy as np
from scipy.constants import R

pluto = {
    'name':                'Pluto',
    'radius':              1189.9e3,
    'mass':                1.3895506e22,
    'constant_obliquity':  0.,
    'surface_temperature': 40.,
    'layers':              {
        'Core':      {
            # Pre-calculated Layer Structure
            'radius_upper':                      910.903e3,
            'radius_lower':                      0.,
            'mass':                              1.0298825e22,
            # Material Parameters
            'material_density':                  3300.,
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
            # Other Model Inputs
            'growth_layer':                      False,
            'constant_ocean_temperature':        0.,
            'constant_viscoelastic_temperature': 0.

        },
        'Icy Shell': {
            # Pre-calculated Layer Structure
            'radius_upper':                      1189.9e3,
            'radius_lower':                      910.903e3,
            'mass':                              3.59668e+21,
            # Material Parameters
            'material_density':                  920.,
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
            # Other Model Inputs
            'growth_layer':                      True,
            'constant_ocean_temperature':        273.15,
            'constant_viscoelastic_temperature': 273.,
            'viscoelastic_top_temperature':      273.15 / ((np.log(10.) / (58.4e3 / (R * 273.15))) + 1.)
        }
    }
}

charon = {
    'name':                'Charon',
    'radius':              606000.0,
    'mass':                1.670299e21,
    'constant_obliquity':  0.,
    'surface_temperature': 40.,
    'layers':              {
        'Core':      {
            # Pre-calculated Layer Structure
            'radius_upper':                      437374.8067,
            'radius_lower':                      0.,
            'mass':                              1.1356611e21,
            # Material Parameters
            'material_density':                  3300.,
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
            # Other Model Inputs
            'growth_layer':                      False,
            'constant_ocean_temperature':        0.,
            'constant_viscoelastic_temperature': 0.

        },
        'Icy Shell': {
            # Pre-calculated Layer Structure
            'radius_upper':                      606000.0,
            'radius_lower':                      437374.8067,
            'mass':                              5.346380e20,
            # Material Parameters
            'material_density':                  920.,
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
            # Other Model Inputs
            'growth_layer':                      True,
            'constant_ocean_temperature':        273.15,
            'constant_viscoelastic_temperature': 273.,
            'viscoelastic_top_temperature':      273.15 / ((np.log(10.) / (58.4e3 / (R * 273.15))) + 1.)
        }
    }
}
