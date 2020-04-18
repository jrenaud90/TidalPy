import numpy as np
from scipy.constants import G
import matplotlib.pyplot as plt

from TidalPy.planets import build_planet

# Computation number
N = 10

# Eventually we can just build the planets with real materials from BurnMan...
# triton = build_planet('Triton')
# neptune = build_planet('Neptune')

# ...But, for testing, lets keep it simple and build the planet inside this script...

# First the host and orbital parameters
host_mass = 1.024e26
host_radius = 24.622e6
eccentricity = 0.01
inclination = 0.0
semi_major_axis = 354759.e3

# Now the planet parameters

radius = 1.3534e6
core_radius = .3 * radius
mantle_radius = .8 * radius
core_density = 5000.
mantle_density = 3000.
crust_density = 1000.
volume = (4. / 3.) * np.pi * radius**3

# Planet spin rate is set by the spin_rate_type and spin_rate_value
#    If the spin_rate_type is set to 'manual' then the spin_rate_value will be used as is, assuming rad/s
#    If the spin_rate_type is set to 'sync' then the spin_rate_value will be ignored and the spin-rate == orbital motion
#    If the spin_rate_type is set to 'frac' then the spin-rate will = spin_rate_value * orbital_motion
spin_rate_type = 'sync'
spin_rate_value = None

# Build the geometry
radius_array_wbottom = np.linspace(0., radius, N)
radius_array = radius_array_wbottom[1:]
num_shells = len(radius_array)
volume_array = (4. / 3.) * np.pi * (radius_array_wbottom[1:]**3 - radius_array_wbottom[:-1]**3)
core_radii = radius_array < core_radius
mantle_radii = np.logical_and(radius_array >= core_radius, radius_array < mantle_radius)
crust_radii = radius_array >= mantle_radius
density_array = np.zeros_like(volume_array)
density_array[core_radii] = core_density
density_array[mantle_radii] = mantle_density
density_array[crust_radii] = crust_density
    
# Determine mass and gravity
mass_array = density_array * volume_array
mass = sum(mass_array)
mass_below_array = np.zeros_like(mass_array)
last_mass_below = 0.
for layer_i, mass in enumerate(mass_array):
    new_mass_below = last_mass_below + mass
    mass_below_array[layer_i] = new_mass_below
    last_mass_below = new_mass_below

gravity_array = G * mass_below_array / radius_array**2
surf_gravity = gravity_array[-1]

# Update orbital parameters now that we have the target planet's mass
from TidalPy.tools.conversions import semi_a2orbital_motion
orbital_motion = semi_a2orbital_motion(semi_major_axis, host_mass, mass)
if spin_rate_type == 'sync':
    spin_rate = orbital_motion
elif spin_rate_type == 'frac':
    spin_rate = spin_rate_value * orbital_motion
elif spin_rate_type == 'manual':
    spin_rate = spin_rate_value
else:
    raise NotImplementedError

# Get a visual on the interior structure
fig, ax_gravity = plt.subplots()
ax_gravity.plot(radius_array/1000, gravity_array, '-g', label='Gravity')
ax_gravity.set_xlabel('Radius [km]')
ax_gravity.set_ylabel('Gravity [m s-2]')

ax_density = ax_gravity.twinx()
ax_density.plot(radius_array/1000, density_array, ':k', label='Density')
ax_density.set_ylabel('Density [kg m-3]')
plt.show()

# Rheological Properties
core_shear = 100.e9
core_visc = 1.e24
mantle_shear = 50.e9
mantle_visc = 1.e20
crust_shear = 3.e9
crust_visc = 1.e18

viscosity_array = np.zeros_like(radius_array)
viscosity_array[core_radii] = core_visc
viscosity_array[mantle_radii] = mantle_visc
viscosity_array[crust_radii] = crust_visc

shear_array = np.zeros_like(radius_array)
shear_array[core_radii] = core_shear
shear_array[mantle_radii] = mantle_shear
shear_array[crust_radii] = crust_shear

# Solve for complex compliance. For now we are only looking at one mode, at the orbital frequency
#    TODO: Add multi-mode calculation
from TidalPy.rheology.complexCompliance.compliance_models import andrade
complex_compliance_array = andrade(orbital_motion, shear_array**(-1), viscosity_array, alpha=0.3, zeta=1.)
complex_shear_array = complex_compliance_array**(-1)

# Multilayer Calculation
from TidalPy.tides.love3d_dev.propagation_matrix import build_fundamental_matrix, find_propagator_matrix

# Build propagator matrix
Y, Y_inv, Y_reduced_shifted = \
    build_fundamental_matrix(radius_array, complex_shear_array, density_array, gravity_array, order_l=2)
# Y_reduced_shifted == Y_i @ Y_i-1 ^ -1

# Build seed matrix at core (central boundary conditions)
seed_matrix = np.zeros((6, 3), dtype=np.complex)
# From IcyDwarf: "Central boundary conditions (3). They are inconsequential on the rest of the solution, so false assumptions are OK."
# seed_matrix[0, 2] = 1.0 + 0.j# Roberts & Nimmo (2008): liquid innermost zone.
# seed_matrix[1, 3] = 1.0 + 0.j
# seed_matrix[2, 5] = 1.0 + 0.j

seed_matrix[0, 0] = 1.0 + 0.j # Alternative: Henning & Hurford (2014): solid innermost zone
seed_matrix[1, 1] = 1.0 + 0.j
seed_matrix[2, 2] = 1.0 + 0.j

# seed_matrix[0, 0] = 0.05 + 0.j # Boundary conditions for Tobie et al. (2005) benchmark
# seed_matrix[1, 1] = 0.01 + 0.j
# seed_matrix[2, 5] = 1.0 + 0.j

# TODO: Confused how this seed matrix is built. from Henning & Hurford
#   "At the core, a special seed matrix Bcore is created with only three columns, equal to the first, second,
#      and third columns of Y for the properties at the base layer"
#   But in icydwarf it looks like it is just set by constant numbers?

# METHOD 2 from H&H2014
seed_matrix_new = np.zeros((6, 3), dtype=np.complex)
seed_matrix_new[:, 0] = Y[:, 0, 0]
seed_matrix_new[:, 1] = Y[:, 1, 0]
seed_matrix_new[:, 2] = Y[:, 2, 0]

# Find the rest of the propagator matrix at all other layers
propagator_matrix = find_propagator_matrix(Y_reduced_shifted, seed_matrix_new)

# Apply surface boundary conditions
# TODO: Note the -5 is actually -(2l + 1), also the negative seems to depend on some other definition (tobie 2005 has positive)
b_surf = np.asarray([0., 0., -5. / radius], dtype=np.complex)

# Surface matrix is the 3rd, 4th, and 6th columns of the aggregate matrix at the top-most shell
surf_matrix = np.zeros((3, 3), dtype=np.complex)
surf_matrix[0, :] = propagator_matrix[2, :, -1]
surf_matrix[1, :] = propagator_matrix[3, :, -1]
surf_matrix[2, :] = propagator_matrix[5, :, -1]

# We need the inverse of this surface matrix
surf_matrix_inv = np.linalg.inv(surf_matrix)

# Now apply the boundary condition and store in vector c
c_vector = surf_matrix_inv @ b_surf

# Calculate the y (tide) vector
y_tide = np.zeros((6, num_shells), dtype=np.complex)
for shell_i in range(num_shells):
    y_tide[:, shell_i] = propagator_matrix[:, :, shell_i] @ c_vector

# Calculate Love and Shida Numbers
love_number = -y_tide[4] + 1 # TODO: Different than Tobie 2005 eq. 9
h_number = gravity_array * y_tide[0]
l_number = gravity_array * y_tide[1]

# Show Love number
fig, ax_love = plt.subplots()
ax_love.plot(radius_array/1000, -np.imag(love_number), '-r')
ax_love.set_xlabel('Radius [km]')
ax_love.set_ylabel('-Im[k_2]')
plt.show()

breakpoint()

# Tidal Potential
#    # TODO

potential = 1
potential_dTheta = 1
potential_dPhi = 1
potential_dTheta2 = 1
potential_dTheta_dPhi = 1
