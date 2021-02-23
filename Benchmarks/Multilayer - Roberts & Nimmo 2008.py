""" Benchmarking TidalPy's Multi-layer code against Roberts + Nimmo 2008

References
----------
RN08 - Roberts and Nimmo (2008; Icarus; DOI: 10.1016/j.icarus.2007.11.010)
"""

import numpy as np
import matplotlib.pyplot as plt
np.seterr(divide='raise')

from TidalPy.constants import G
from TidalPy.tools.conversions import orbital_motion2semi_a
from TidalPy.utilities.numpy_helper import find_nearest
from TidalPy.rheology.complex_compliance.compliance_models import maxwell_array, sundberg_array
# Load TidalPy's multilayer functions
from TidalPy.tides.multilayer import fundamental_matrix_orderl2, propagate, decompose
from TidalPy.tides.potential.synchronous_low_e import tidal_potential
from TidalPy.tides.multilayer.heating import calc_radial_tidal_heating
from TidalPy.tides.multilayer.stress_strain import calculate_strain


# Planet Structure
R_planet = 1600.e3
rho_mantle = 3300.
# Two different core density models, lower value assumes Fe / Fe-S core composition.
rho_core_m1 = 5150.
rho_core_m2 = 8000.
rho_average = 3500.

viscosity_mantle = 1.e20
shear_mantle = 67.e9

# From RN08 "Rigidity and viscosity have been reduced by a factor of 1e6 and 1e9 respectively from the overlying ice shell"
viscosity_core = 1.e13 / 1.e9
shear_core = 4.e9 / 1.e6

# Core radius is determined from the core density and planet radius
R_core_m1 = ((rho_average - rho_mantle) / (rho_core_m1 - rho_mantle))**(1 / 3) * R_planet
R_core_m2 = ((rho_average - rho_mantle) / (rho_core_m2 - rho_mantle))**(1 / 3) * R_planet

# Orbital properties
planet_mass = 1.08e20
host_mass = 5.683e26
eccentricity = 0.0045
orbital_freq_HZ = 5.308e-5  # This is reported in Hz in RN08. We will need to convert them to rad s-1
orbital_freq = 2. * np.pi * orbital_freq_HZ
semi_major_axis = orbital_motion2semi_a(orbital_freq, host_mass, planet_mass)

# Setup homogeneous domain
N = 500
radius_array = np.linspace(0., R_planet, N+1)[1:]
volume_array = np.zeros_like(radius_array)
volume_array[0] = (4. / 3.) * np.pi * radius_array[0]**3
volume_array[1:] = (4. / 3.) * np.pi * (radius_array[1:]**3 - radius_array[:-1]**3)

# There are three different models we are comparing:
#    m0 = Homogenious interior with average density.
#    m1 = Differentiated interior with a liquid core at 5150 kg m-3 and a mantle at 3300 kg m-3
#    m2 = Differentiated interior with a liquid core at 8000 kg m-3 and a mantle at 3300 kg m-3
density_array_m0 = rho_average * np.ones_like(radius_array)
density_array_m1 = rho_mantle * np.ones_like(radius_array)
density_array_m1[radius_array <= R_core_m1] = rho_core_m1
density_array_m2 = rho_mantle * np.ones_like(radius_array)
density_array_m2[radius_array <= R_core_m2] = rho_core_m2

# Calculate masses, mass below, gravities
mass_array_m0 = volume_array * density_array_m0
mass_array_m1 = volume_array * density_array_m1
mass_array_m2 = volume_array * density_array_m2

masses_below_m0 = np.zeros_like(radius_array)
masses_below_m0[0] = mass_array_m0[0]
masses_below_m0[1:] = np.asarray([sum(mass_array_m0[0:i+1]) for i in range(1, N)])
masses_below_m1 = np.zeros_like(radius_array)
masses_below_m1[0] = mass_array_m1[0]
masses_below_m1[1:] = np.asarray([sum(mass_array_m1[0:i+1]) for i in range(1, N)])
masses_below_m2 = np.zeros_like(radius_array)
masses_below_m2[0] = mass_array_m2[0]
masses_below_m2[1:] = np.asarray([sum(mass_array_m2[0:i+1]) for i in range(1, N)])

gravity_array_m0 = G * masses_below_m0 / radius_array**2
gravity_array_m1 = G * masses_below_m1 / radius_array**2
gravity_array_m2 = G * masses_below_m2 / radius_array**2

gravity_array_by_model = {'HG': gravity_array_m0, 'LC0': gravity_array_m1, 'LC1': gravity_array_m2}
mass_array_by_model = {'HG': mass_array_m0, 'LC0': mass_array_m1, 'LC1': mass_array_m2}
density_array_by_model = {'HG': density_array_m0, 'LC0': density_array_m1, 'LC1': density_array_m2}

# Plot interior models
line_style_by_model = {'HG': '-', 'LC0': '--', 'LC1': ':'}
fig_int, axes_int = plt.subplots(ncols=3)
ax_int_rho = axes_int[0]
ax_int_mass = axes_int[1]
ax_int_grav = axes_int[2]

ax_int_rho.set(ylabel='Radius [km]', xlabel='Density [kg m-3]')
ax_int_mass.set(xlabel='Mass [kg]')
ax_int_grav.set(xlabel='Gravity Acc. [m s-2]')
fig_int.suptitle('Interior Models')

for model_name, model_ls in line_style_by_model.items():
    ax_int_rho.plot(density_array_by_model[model_name], radius_array/1000., label=model_name, c='k', ls=model_ls)
    ax_int_mass.plot(mass_array_by_model[model_name], radius_array/1000., label=model_name, c='k', ls=model_ls)
    ax_int_grav.plot(gravity_array_by_model[model_name], radius_array/1000., label=model_name, c='k', ls=model_ls)

ax_int_rho.legend(loc='best')
fig_int.tight_layout()
plt.show()

# Calculate rheological properties
viscosity_array_m0 = viscosity_mantle * np.ones_like(radius_array)
viscosity_array_m1 = viscosity_mantle * np.ones_like(radius_array)
viscosity_array_m1[radius_array <= R_core_m1] = viscosity_core
viscosity_array_m2 = viscosity_mantle * np.ones_like(radius_array)
viscosity_array_m2[radius_array <= R_core_m2] = viscosity_core

shear_array_m0 = shear_mantle * np.ones_like(radius_array)
shear_array_m1 = shear_mantle * np.ones_like(radius_array)
shear_array_m1[radius_array <= R_core_m1] = shear_core
shear_array_m2 = shear_mantle * np.ones_like(radius_array)
shear_array_m2[radius_array <= R_core_m2] = shear_core

# Use the Maxwell rheology
complex_compliance_array_m0 = maxwell_array(orbital_freq, shear_array_m0**(-1), viscosity_array_m0)
complex_compliance_array_m1 = maxwell_array(orbital_freq, shear_array_m1**(-1), viscosity_array_m1)
complex_compliance_array_m2 = maxwell_array(orbital_freq, shear_array_m2**(-1), viscosity_array_m2)
complex_shear_array_m0 = complex_compliance_array_m0**(-1)
complex_shear_array_m1 = complex_compliance_array_m1**(-1)
complex_shear_array_m2 = complex_compliance_array_m2**(-1)

# Find the fundamental matrix; skip the innermost shell as that will be a boundary condition.
Y_m0, Y_inv_m0, derivative_mtx_m0 = \
    fundamental_matrix_orderl2(radius_array, complex_shear_array_m0, density_array_m0, gravity_array_m0)
Y_m1, Y_inv_m1, derivative_mtx_m1 = \
    fundamental_matrix_orderl2(radius_array[radius_array > R_core_m1],
                               complex_shear_array_m1[radius_array > R_core_m1],
                               density_array_m1[radius_array > R_core_m1],
                               gravity_array_m1[radius_array > R_core_m1])
Y_m2, Y_inv_m2, derivative_mtx_m2 = \
    fundamental_matrix_orderl2(radius_array[radius_array > R_core_m2],
                               complex_shear_array_m2[radius_array > R_core_m2],
                               density_array_m2[radius_array > R_core_m2],
                               gravity_array_m2[radius_array > R_core_m2])

# Propagate the tidal solution through the world's shells
central_boundary_condition = np.zeros((6, 3))
# Roberts & Nimmo (2008): Liquid innermost zone.
central_boundary_condition[0, 0] = 0.05
central_boundary_condition[1, 1] = 0.01
central_boundary_condition[5, 2] = 1.

tidal_y_m0, tidal_y_deriv_m0 = \
    propagate(Y_m0, Y_inv_m0, derivative_mtx_m0, central_boundary_condition, R_planet, order_l=2)
tidal_y_m1, tidal_y_deriv_m1 = \
    propagate(Y_m1, Y_inv_m1, derivative_mtx_m1, central_boundary_condition, R_planet, order_l=2)
tidal_y_m2, tidal_y_deriv_m2 = \
    propagate(Y_m2, Y_inv_m2, derivative_mtx_m2, central_boundary_condition, R_planet, order_l=2)
tidal_y_by_model = {'HG': tidal_y_m0, 'LC0': tidal_y_m1, 'LC1': tidal_y_m2}
solved_radii_by_model = {'HG': radius_array, 'LC0': radius_array[radius_array > R_core_m1], 'LC1': radius_array[radius_array > R_core_m2]}


# Plot tidal y's
fig_tidal_y, axes_tidal_y = plt.subplots(ncols=4, figsize=(10, 5))
ax_y1 = axes_tidal_y[0]
ax_y3 = axes_tidal_y[1]
ax_y2 = axes_tidal_y[2]
ax_y4 = axes_tidal_y[3]

ax_y1.set(ylabel='Radius [km]', xlabel='$y_{1}$ [m / (m/s)$^{2}$]', xlim=(0, 0.15))
ax_y3.set(ylabel='Radius [km]', xlabel='$y_{3}$ [m / (m/s)$^{2}$]', xlim=(-0.04, 0.04))
ax_y2.set(ylabel='Radius [km]', xlabel='$y_{2}$ [kg / m$^{3}$]', xlim=(-2000, 5000))
ax_y4.set(ylabel='Radius [km]', xlabel='$y_{4}$ [kg / m$^{3}$]', xlim=(0, 2000))

# Load in RN08 data
RN08_C = 'b'
T05_C = 'r'
plot_RN08 = True
plot_T05 = True

if plot_RN08:
    RN08_data = np.loadtxt('RN08-Data.csv', skiprows=1, delimiter=',', dtype=str)
    RN08_data[RN08_data == ''] = np.nan
    RN08_data = RN08_data.astype(np.float)
    ax_y1.scatter(RN08_data[:, 0], RN08_data[:, 1]/1000., label='RN08-HG', c=RN08_C, marker='.', s=20)
    ax_y1.scatter(RN08_data[:, 2], RN08_data[:, 3]/1000., label='RN08-LC0', c=RN08_C, marker='+', s=20)

    ax_y3.scatter(RN08_data[:, 4], RN08_data[:, 5]/1000., label='RN08-HG', c=RN08_C, marker='.', s=20)
    ax_y3.scatter(RN08_data[:, 6], RN08_data[:, 7]/1000., label='RN08-LC0', c=RN08_C, marker='+', s=20)

    ax_y2.scatter(RN08_data[:, 8], RN08_data[:, 9]/1000., label='RN08-HG', c=RN08_C, marker='.', s=20)
    ax_y2.scatter(RN08_data[:, 10], RN08_data[:, 11]/1000., label='RN08-LC0', c=RN08_C, marker='+', s=20)

    ax_y4.scatter(RN08_data[:, 12], RN08_data[:, 13]/1000., label='RN08-HG', c=RN08_C, marker='.', s=20)
    ax_y4.scatter(RN08_data[:, 14], RN08_data[:, 15]/1000., label='RN08-LC0', c=RN08_C, marker='+', s=20)

if plot_T05:
    T05_data = np.loadtxt('T05-Data.csv', skiprows=1, delimiter=',', dtype=str)
    T05_data[T05_data == ''] = np.nan
    T05_data = T05_data.astype(np.float)
    ax_y1.scatter(T05_data[:, 0], T05_data[:, 1] / 1000., label='T05-HG', c=T05_C, marker='.', s=20)
    ax_y1.scatter(T05_data[:, 2], T05_data[:, 3] / 1000., label='T05-LC0', c=T05_C, marker='+', s=20)
    ax_y1.scatter(T05_data[:, 4], T05_data[:, 5] / 1000., label='T05-LC1', c=T05_C, marker='1', s=20)

    ax_y3.scatter(T05_data[:, 6], T05_data[:, 7] / 1000., label='T05-HG', c=T05_C, marker='.', s=20)
    ax_y3.scatter(T05_data[:, 8], T05_data[:, 9] / 1000., label='T05-LC0', c=T05_C, marker='+', s=20)
    ax_y3.scatter(T05_data[:, 10], T05_data[:, 11] / 1000., label='T05-LC1', c=T05_C, marker='1', s=20)

    ax_y2.scatter(T05_data[:, 12], T05_data[:, 13] / 1000., label='T05-HG', c=T05_C, marker='.', s=20)
    ax_y2.scatter(T05_data[:, 14], T05_data[:, 15] / 1000., label='T05-LC0', c=T05_C, marker='+', s=20)
    ax_y2.scatter(T05_data[:, 16], T05_data[:, 17] / 1000., label='T05-LC1', c=T05_C, marker='1', s=20)

    ax_y4.scatter(T05_data[:, 18], T05_data[:, 19] / 1000., label='T05-HG', c=T05_C, marker='.', s=20)
    ax_y4.scatter(T05_data[:, 20], T05_data[:, 21] / 1000., label='T05-LC0', c=T05_C, marker='+', s=20)
    ax_y4.scatter(T05_data[:, 22], T05_data[:, 23] / 1000., label='T05-LC1', c=T05_C, marker='1', s=20)

# Plot our results
for model_name, model_ls in line_style_by_model.items():
    tidal_y_set = tidal_y_by_model[model_name]
    radii = solved_radii_by_model[model_name]
    ax_y1.plot((tidal_y_set[0, :]), radii/1000., label=model_name, c='k', ls=model_ls)
    ax_y3.plot((tidal_y_set[2, :]), radii/1000., label=model_name, c='k', ls=model_ls)
    ax_y2.plot((tidal_y_set[1, :]), radii/1000., label=model_name, c='k', ls=model_ls)
    ax_y4.plot((tidal_y_set[3, :]), radii/1000., label=model_name, c='k', ls=model_ls)

ax_y1.legend(loc='best')
fig_tidal_y.tight_layout()
plt.show()