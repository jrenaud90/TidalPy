import numpy as np
import matplotlib.pyplot as plt


from TidalPy.constants import G
from TidalPy.tools.conversions import orbital_motion2semi_a
from TidalPy.utilities.numpy_helper import find_nearest
from TidalPy.rheology.complex_compliance.compliance_models import maxwell_array, sundberg_array
# Load TidalPy's multilayer functions
from TidalPy.tides.multilayer import fundamental_matrix_orderl2, propagate, decompose
from TidalPy.tides.potential.synchronous_low_e import tidal_potential
from TidalPy.tides.multilayer.heating import calc_radial_tidal_heating
from TidalPy.tides.multilayer.stress_strain import calculate_strain
from tides.multilayer.numerical_int.radial_functions_solid import radial_func_derivatives_general

full_domain = True

# Planet Structure
R_planet = 1600.e3
rho_mantle = 3300.
# Two different core density models, lower value assumes Fe / Fe-S core composition.
rho_core_m1 = 5150.
rho_core_m2 = 8000.
rho_average = 3500.

use_incompressibility = False # To compare to RN08

viscosity_mantle = 1.e20
shear_mantle_0 = 6.68250e10
shear_mantle_1 = 6.68250e10
shear_mantle_2 = 6.68250e10
if use_incompressibility:
    bulk_mantle_0 = 1.e50
    bulk_mantle_1 = 1.e50
    bulk_mantle_2 = 1.e50
else:
    bulk_mantle_0 = 1.2210e11
    bulk_mantle_1 = 1.2210e11
    bulk_mantle_2 = 1.2210e11

# From RN08 "Rigidity and viscosity have been reduced by a factor of 1e6 and 1e9 respectively from the overlying ice shell"
viscosity_core = 1.e13 / 1.e9
shear_core_0 = 0.
shear_core_1 = 6.68250e-5
shear_core_2 = 6.68250e-5

if use_incompressibility:
    bulk_core_0 = 0.
    bulk_core_1 = 1.e50
    bulk_core_2 = 1.e50
else:
    bulk_core_0 = 0.
    bulk_core_1 = 1.8e11
    bulk_core_2 = 2.88e11

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
N = 100
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

# Calculate rheological properties
viscosity_array_m0 = viscosity_mantle * np.ones_like(radius_array)
viscosity_array_m1 = viscosity_mantle * np.ones_like(radius_array)
viscosity_array_m1[radius_array <= R_core_m1] = viscosity_core
viscosity_array_m2 = viscosity_mantle * np.ones_like(radius_array)
viscosity_array_m2[radius_array <= R_core_m2] = viscosity_core

shear_array_m0 = shear_mantle_0 * np.ones_like(radius_array)
shear_array_m1 = shear_mantle_1 * np.ones_like(radius_array)
shear_array_m1[radius_array <= R_core_m1] = shear_core_1
shear_array_m2 = shear_mantle_2 * np.ones_like(radius_array)
shear_array_m2[radius_array <= R_core_m2] = shear_core_2

bulk_array_m0 = bulk_mantle_0 * np.ones_like(radius_array)
bulk_array_m1 = bulk_mantle_1 * np.ones_like(radius_array)
bulk_array_m1[radius_array <= R_core_m1] = bulk_core_1
bulk_array_m2 = bulk_mantle_2 * np.ones_like(radius_array)
bulk_array_m2[radius_array <= R_core_m2] = bulk_core_2

# Use the Maxwell rheology
complex_compliance_array_m0 = maxwell_array(orbital_freq, shear_array_m0**(-1), viscosity_array_m0)
complex_compliance_array_m1 = maxwell_array(orbital_freq, shear_array_m1**(-1), viscosity_array_m1)
complex_compliance_array_m2 = maxwell_array(orbital_freq, shear_array_m2**(-1), viscosity_array_m2)
complex_shear_array_m0 = complex_compliance_array_m0**(-1)
complex_shear_array_m1 = complex_compliance_array_m1**(-1)
complex_shear_array_m2 = complex_compliance_array_m2**(-1)



interp_funcs = {
    'm0': (lambda x: np.interp(x, radius_array, complex_shear_array_m0),
           lambda x: np.interp(x, radius_array, bulk_array_m0),
           lambda x: np.interp(x, radius_array, density_array_m0),
           lambda x: np.interp(x, radius_array, gravity_array_m0)),
    'm1': (lambda x: np.interp(x, radius_array, complex_shear_array_m1),
           lambda x: np.interp(x, radius_array, bulk_array_m1),
           lambda x: np.interp(x, radius_array, density_array_m1),
           lambda x: np.interp(x, radius_array, gravity_array_m1)),
    'm2': (lambda x: np.interp(x, radius_array, complex_shear_array_m2),
           lambda x: np.interp(x, radius_array, bulk_array_m2),
           lambda x: np.interp(x, radius_array, density_array_m2),
           lambda x: np.interp(x, radius_array, gravity_array_m2))
}

def deriv_m0(x, y):
    shear_modulus, bulk_modulus, density, gravity = interp_funcs['m0']

    result = radial_func_derivatives_general(x, y, shear_modulus(x), bulk_modulus(x), density(x), gravity(x), orbital_freq)

    return result

def deriv_m1(x, y):
    shear_modulus, bulk_modulus, density, gravity = interp_funcs['m1']

    result = radial_func_derivatives_general(x, y, shear_modulus(x), bulk_modulus(x), density(x), gravity(x), orbital_freq)

    return result

def deriv_m2(x, y):
    shear_modulus, bulk_modulus, density, gravity = interp_funcs['m2']

    result = radial_func_derivatives_general(x, y, shear_modulus(x), bulk_modulus(x), density(x), gravity(x), orbital_freq)

    return result

surf_condition = (5. / R_planet)

def bc_residual(ya, yb):

    # At the core
    # y1 = y3 = y5 = 0

    # At the surface
    # y2, y4, y6 = (0, 0, 2l + 1 / R)

    return (ya[0], yb[1], ya[2], yb[3], ya[4], yb[5] - surf_condition)


from scipy.integrate import solve_bvp

if full_domain:
    solved_radii_by_model = {'HG' : radius_array, 'LC0': radius_array,
                             'LC1': radius_array}
else:
    solved_radii_by_model = {'HG' : radius_array, 'LC0': radius_array[radius_array > R_core_m1],
                             'LC1': radius_array[radius_array > R_core_m2]}

# Initial Guesses...
initial_guesses = {}
for structure_name, radii in solved_radii_by_model.items():

    y1 = np.ones_like(radii)
    y2 = np.ones_like(radii)
    y3 = np.ones_like(radii)
    y4 = np.ones_like(radii)
    y5 = np.ones_like(radii)
    y6 = np.ones_like(radii)

    guess = np.vstack((y1, y2, y3, y4, y5, y6))
    initial_guesses[structure_name] = guess

tidal_y_by_model = dict()
solved_radii_by_model_bvp = dict()
for structure_name, deriv_func in (('HG', deriv_m0), ('LC0', deriv_m1), ('LC1', deriv_m2)):

    print(f'Working on {structure_name}.')

    radii = solved_radii_by_model[structure_name]
    initial_guess = initial_guesses[structure_name]
    solution = solve_bvp(deriv_func, bc_residual, radii, initial_guess, tol=1e-3, verbose=2)

    if solution.status != 0:
        print("WARNING: sol.status is %d" % solution.status)
    print(solution.message)

    tidal_y_by_model[structure_name] = solution.y
    solved_radii_by_model_bvp[structure_name] = solution.x



# Plot tidal y's
fig_tidal_y, axes_tidal_y = plt.subplots(ncols=4, figsize=(10, 5))
ax_y1 = axes_tidal_y[0]
ax_y3 = axes_tidal_y[1]
ax_y2 = axes_tidal_y[2]
ax_y4 = axes_tidal_y[3]
line_style_by_model = {'HG': '-', 'LC0': '--', 'LC1': ':'}

ax_y1.set(ylabel='Radius [km]', xlabel='$y_{1}$ [m / (m/s)$^{2}$]', xlim=(0, 0.15))
ax_y3.set(ylabel='Radius [km]', xlabel='$y_{3}$ [m / (m/s)$^{2}$]', xlim=(-0.04, 0.04))
ax_y2.set(ylabel='Radius [km]', xlabel='$y_{2}$ [kg / m$^{3}$]', xlim=(-2000, 8000))
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
    radii = solved_radii_by_model_bvp[model_name]
    ax_y1.plot((tidal_y_set[0, :]), radii/1000., label=model_name, c='k', ls=model_ls)
    ax_y3.plot((tidal_y_set[2, :]), radii/1000., label=model_name, c='k', ls=model_ls)
    ax_y2.plot((tidal_y_set[1, :]), radii/1000., label=model_name, c='k', ls=model_ls)
    ax_y4.plot((tidal_y_set[3, :]), radii/1000., label=model_name, c='k', ls=model_ls)

ax_y1.legend(loc='best')
fig_tidal_y.tight_layout()
plt.show()