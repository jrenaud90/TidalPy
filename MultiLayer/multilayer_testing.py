import numpy as np
from scipy.constants import G
import matplotlib.pyplot as plt

from TidalPy import build_world

# Build Io, but leave the host just as numbers for now.
io = build_world('Earth')
io.paint()
host_mass = 1.89813e27
semi_major_axis = 421700.0e3

breakpoint()

# First the host and orbital parameters
host_mass = 1.024e26
host_radius = 24.622e6
eccentricity = 0.1
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
N = 100
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

# Get a visual on the interior structure
fig, ax_gravity = plt.subplots()
ax_gravity.plot(radius_array/1000, gravity_array, '-g', label='Gravity')
ax_gravity.set_xlabel('Radius [km]')
ax_gravity.set_ylabel('Gravity [m s-2]')

ax_density = ax_gravity.twinx()
ax_density.plot(radius_array/1000, density_array, ':k', label='Density')
ax_density.set_ylabel('Density [kg m-3]')
plt.show()

# We do not want the inner-most shell that TidalPy builds as the r=0 is taken into account in the boundary conditions.
# radii = io.radii[1:]
# gravities = io.gravities[1:]
# densities = io.densities[1:]
# num_shells = len(io.radii[1:])
# radius = io.radius
# mass = io.mass
# surf_gravity = gravities[-1]
# core_radii = radii < io.core.radius
# mantle_radii = radii >= io.core.radius

radii = radius_array
gravities = gravity_array
densities = density_array
num_shells = num_shells
radius = radius

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
forcing_frequency = orbital_motion

# Rheological Properties - Hardcoded for now
core_bulk = 150.e9
core_shear = 100.e9
core_visc = 1.e24

mantle_bulk = 131.e9
mantle_shear = 50.e9
mantle_visc = 1.e20

crust_bulk = 125.e9
crust_shear = 50.e9
crust_visc = 1.e19

viscosity_array = np.zeros_like(radii)
viscosity_array[core_radii] = core_visc
viscosity_array[mantle_radii] = mantle_visc
viscosity_array[crust_radii] = crust_visc

shear_array = np.zeros_like(radii)
shear_array[core_radii] = core_shear
shear_array[mantle_radii] = mantle_shear
shear_array[crust_radii] = crust_shear

bulk_array = np.zeros_like(radii)
bulk_array[core_radii] = core_bulk
bulk_array[mantle_radii] = mantle_bulk
bulk_array[crust_radii] = crust_bulk

# Solve for complex compliance. For now we are only looking at one mode, at the orbital frequency
from TidalPy.rheology.complex_compliance.compliance_models import maxwell_array

complex_compliance_array = maxwell_array(forcing_frequency, shear_array**(-1), viscosity_array)
complex_shear_array = complex_compliance_array**(-1)

# Multilayer Calculation
from TidalPy.tides.love3d_dev.propagation_matrix import build_fundamental_matrix, calculate_compound_matrix

# Build propagator matrix
Y, Y_inv, Y_reduced_shifted = \
    build_fundamental_matrix(radii, complex_shear_array, densities, gravities, order_l=2)
# Y_reduced_shifted == Y_i @ Y_i-1 ^ -1

# Build seed matrix at core (central boundary conditions)
seed_matrix_other = np.zeros((6, 6), dtype=np.complex)
# From IcyDwarf: "Central boundary conditions (3). They are inconsequential on the rest of the solution, so false assumptions are OK."
seed_matrix_other[2, 0] = 1.0 + 0.j # Roberts & Nimmo (2008): liquid innermost zone.
seed_matrix_other[3, 1] = 1.0 + 0.j
seed_matrix_other[5, 2] = 1.0 + 0.j

# seed_matrix_other[0, 0] = 1.0 + 0.j  # Alternative: Henning & Hurford (2014): solid innermost zone
# seed_matrix_other[1, 1] = 1.0 + 0.j
# seed_matrix_other[2, 2] = 1.0 + 0.j

# seed_matrix_other[0, 0] = 0.05 + 0.j # Boundary conditions for Tobie et al. (2005) benchmark
# seed_matrix_other[1, 1] = 0.01 + 0.j
# seed_matrix_other[5, 2] = 1.0 + 0.j

# TODO: Confused how this seed matrix is built. from Henning & Hurford
#   "At the core, a special seed matrix Bcore is created with only three columns, equal to the first, second,
#      and third columns of Y for the properties at the base layer"
#   But in icydwarf it looks like it is just set by constant numbers?

# METHOD 2 from H&H2014
# seed_matrix_new = np.zeros((6, 3), dtype=np.complex)
# seed_matrix_new[:, 0] = Y[:, 0, 0]
# seed_matrix_new[:, 1] = Y[:, 1, 0]
# seed_matrix_new[:, 2] = Y[:, 2, 0]

# Find the rest of the propagator matrix at all other layers
seed_matrix = seed_matrix_other
compound_matrix = calculate_compound_matrix(Y, Y_inv, seed_matrix)

# Apply surface boundary conditions
# Surface matrix is the 3rd, 4th, and 6th columns of the aggregate matrix at the top-most shell
# TODO: But IcyDwarf defines it as "Define Mbc = 3x3 matrix, rows 3, 4, 6 of Bpropmtx[NR-1]"
surf_matrix = np.zeros((3, 3), dtype=np.complex)
surf_matrix[0, 0] = compound_matrix[2, 0, -1]
surf_matrix[1, 0] = compound_matrix[3, 0, -1]
surf_matrix[2, 0] = compound_matrix[5, 0, -1]

surf_matrix[0, 1] = compound_matrix[2, 1, -1]
surf_matrix[1, 1] = compound_matrix[3, 1, -1]
surf_matrix[2, 1] = compound_matrix[5, 1, -1]

surf_matrix[0, 2] = compound_matrix[2, 2, -1]
surf_matrix[1, 2] = compound_matrix[3, 2, -1]
surf_matrix[2, 2] = compound_matrix[5, 2, -1]

# TODO: Note the -5 is actually -(2l + 1), also the negative seems to depend on some other definition (tobie 2005 has positive)
b_surf = np.asarray([0., 0., -5. / radius], dtype=np.complex)

# We need the inverse of this surface matrix
surf_matrix_inv = np.linalg.inv(surf_matrix)

# Now apply the boundary condition and store in vector c
c_vector = surf_matrix_inv @ b_surf

# Calculate the y (tide) vector
y_tide = np.zeros((6, num_shells), dtype=np.complex)
for shell_i in range(num_shells):
    y_tide[:, shell_i] = compound_matrix[:, :3, shell_i] @ c_vector

# Calculate Love and Shida Numbers
k_number = -y_tide[4] + 1
# k2 different between Tobie2005, TODO: Different than Tobie 2005 eq. 9
h_number = gravities * y_tide[0]
l_number = gravities * y_tide[1]

# Show Love number
fig0, ax0 = plt.subplots()
ax0.plot(radii / 1000, y_tide[0, :], '-r')
ax0.set_xlabel('Radius [km]')
ax0.set_ylabel('Radial Displacement [m-1 s2]')
plt.show()

fig1, ax1 = plt.subplots()
ax1.plot(radii / 1000, y_tide[1, :], '-r')
ax1.set_xlabel('Radius [km]')
ax1.set_ylabel('Tangential Displacement [m-1 s2]')
plt.show()

fig2, ax2 = plt.subplots()
ax2.plot(radii / 1000, y_tide[2, :], '-r')
ax2.set_xlabel('Radius [km]')
ax2.set_ylabel('Radial Stress [kg m-3]')
plt.show()

fig3, ax3 = plt.subplots()
ax3.plot(radii / 1000, y_tide[3, :], '-r')
ax3.set_xlabel('Radius [km]')
ax3.set_ylabel('Tangential Stress [kg m-3]')
plt.show()

fig4, ax4 = plt.subplots()
ax4.plot(radii / 1000, -np.imag(k_number), '-r')
ax4.set_xlabel('Radius [km]')
ax4.set_ylabel('-Im[k_2]')
plt.show()

# Tobie-2005 Heating Rate
y_0_conj_derivative = (np.conjugate(y_tide[0, 1:]) - np.conjugate(y_tide[0, :-1])) / (radii[1:] - radii[:-1])
order_l = 2
sensativity_to_shear = \
    (4. / 3.) * (radii[1:]**2 / np.abs(bulk_array[1:] + (4. / 3.) * complex_shear_array[1:])**2) * np.abs(
            y_tide[2, 1:] - (((bulk_array[1:] - (2. / 3.) * complex_shear_array[1:]) / radii[1:]) *
                            (2. * y_tide[0, 1:] - order_l * (order_l + 1) * y_tide[1, 1:])))**2 - \
    (4. / 3.) * radii[1:] * np.real(y_0_conj_derivative *
                                    (2. * y_tide[0, 1:] - order_l * (order_l + 1) * y_tide[1, 1:])) + \
    (1. / 3.) * np.abs(2. * y_tide[0, 1:] - order_l * (order_l + 1) * y_tide[1, 1:])**2 + \
    order_l * (order_l + 1) * radii[1:]**2 * np.abs(y_tide[3, 1:])**2 / np.abs(complex_shear_array[1:])**2 + \
    order_l * (order_l**2 - 1) * (order_l + 2) * np.abs(y_tide[1, 1:])**2

tidal_heating_per_shell = \
    (21. / 10.) * (forcing_frequency**5 * radius**4 * eccentricity**2 / radii[1:]**2) * \
        sensativity_to_shear * np.imag(complex_shear_array[1:])

fig5, ax5 = plt.subplots()
ax5.plot(radii[1:] / 1000, tidal_heating_per_shell, '-r')
ax5.set_xlabel('Radius [km]')
ax5.set_ylabel('Volumetric Tidal Heating [W m-3]')
plt.show()

print(f'Total Tidal Heating {sum(tidal_heating_per_shell * volume_array[1:]) / 1e12} [TW]')

# Tidal Potential
#    # TODO

potential = 1
potential_dTheta = 1
potential_dPhi = 1
potential_dTheta2 = 1
potential_dTheta_dPhi = 1