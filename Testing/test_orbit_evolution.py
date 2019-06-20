import numpy as np
import matplotlib.pyplot as plt

import TidalPy
from TidalPy.utilities.conversions import days2rads, orbital_motion2semi_a, sec2myr
from TidalPy.utilities.numpy_help import neg_array_for_log_plot

planet = TidalPy.build_planet('charon_example', force_build=True)
host = TidalPy.build_planet('pluto_example', force_build=True)
star = TidalPy.build_planet('sol', force_build=True)

planet.paint(auto_show=True)
host.paint(auto_show=True)

# Inputs
modern_semi_major_axis = np.asarray([1.9591e7])
eccen = np.asarray([0.2])
inclin = np.asarray([5.0])
spin_frequency = np.asarray(days2rads(6.3872304))

# Setup Material Properties
planet.mantle.temperature = 250.
planet.core.set_strength(viscosity=1.e20, shear_modulus=5.e10)
host.mantle.temperature = 250.
host.core.set_strength(viscosity=1.e20, shear_modulus=5.e10)

# Domain of study
semi_major_axis_scale = np.linspace(.2, 2.5, 1000)
# Make sure that the semi-major axis actually has the resonances
resonances = np.asarray(
        [spin_frequency, 2. * spin_frequency, spin_frequency / 2., (2 / 3) * spin_frequency, (3 / 2) * spin_frequency])
resonances = orbital_motion2semi_a(resonances, host.mass, planet.mass) / modern_semi_major_axis
semi_major_axis_scale = np.concatenate((semi_major_axis_scale, resonances))
semi_major_axis_scale = np.sort(semi_major_axis_scale)
semi_major_axis = semi_major_axis_scale * modern_semi_major_axis

# Initialize Orbit & Spin
planet.spin_freq = spin_frequency.copy() * np.ones_like(semi_major_axis)
host.spin_freq = spin_frequency.copy() * np.ones_like(semi_major_axis)
orbit = TidalPy.Orbit(star, host, [planet], duel_dissipation=True, time_study=True)
orbit.set_orbit(planet, semi_major_axis=semi_major_axis, eccentricity=eccen, inclination=inclin,
                inclination_in_deg=True)

### Plot Set 1
print('Plot Set 1 - Individual Evolution')
print('Light colors in the plots indicate negative values.')
for obj in [host, planet]:

    tidal_heating = obj.tidal_heating
    tidal_torque = obj.tidal_ztorque
    spin_derivative = obj.derivative_spin
    tau_spin = -np.abs(obj.spin_freq - orbit.get_orbital_freq(obj)) / spin_derivative

    # Make two tidal torques for log-scale plotting of negative numbers
    tidal_torque_pos, tidal_torque_neg = neg_array_for_log_plot(tidal_torque)

    # Likewise for tau_spin
    tau_spin_pos, tau_spin_neg = neg_array_for_log_plot(tau_spin)

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    axes[0].set_title(obj.name.title())
    axes[0].set_xlabel('Semi-Major Axis (Modern Scale)')
    axes[1].set_xlabel('Semi-Major Axis (Modern Scale)')
    axes[2].set_xlabel('Semi-Major Axis (Modern Scale)')
    axes[0].set_yscale('log')
    axes[1].set_yscale('log')
    axes[2].set_yscale('log')
    axes[0].set_ylabel('Tidal Heating [N m]')
    axes[1].set_ylabel('Tidal Torque (abs) [N m]')
    axes[2].set_ylabel('Spin-Sync Timescale (abs) [Myr]')

    axes[0].plot(semi_major_axis_scale, tidal_heating, c='k')

    axes[1].plot(semi_major_axis_scale, tidal_torque_pos, c='k')
    axes[1].plot(semi_major_axis_scale, tidal_torque_neg, c='k', alpha=0.3)

    axes[2].plot(semi_major_axis_scale, tau_spin_pos, c='k')
    axes[2].plot(semi_major_axis_scale, tau_spin_neg, c='k', alpha=0.3)

plt.show()

# Plot Set 2
print('Plot Set 2 - Orbital Evolution')
print('Light colors in the plots indicate negative values.')

# The evolution of eccentricity and semi-major axis is the same for both the host and the target planet. So we can use
#    either as the reference.
ref_obj = host
de_dt = orbit.get_derivative_eccentricity(ref_obj)
da_dt = orbit.get_derivative_semi_major(ref_obj)

# Get timescales for de_dt
tau_e = -(eccen / de_dt)
tau_e = sec2myr(tau_e)

# Convert units on semi-major axis [m s-1] to [km myr-1]
da_dt = da_dt * (3.154e13 / 1000.)

# Make both positive for log plotting
tau_e_pos, tau_e_neg = neg_array_for_log_plot(tau_e)
da_dt_pos, da_dt_neg = neg_array_for_log_plot(da_dt)

fig2, axes2 = plt.subplots(1, 2, figsize=(18, 9))
axes2[0].set_xlabel('Semi-Major Axis (Modern Scale)')
axes2[1].set_xlabel('Semi-Major Axis (Modern Scale)')
axes2[0].set_yscale('log')
axes2[1].set_yscale('log')
axes2[0].set_ylabel('Circularization Timescale [Myr]')
axes2[1].set_ylabel('Change in Semi-major Axis [km Myr-1]')

axes2[0].plot(semi_major_axis_scale, tau_e_pos, c='k')
axes2[0].plot(semi_major_axis_scale, tau_e_neg, c='k', alpha=0.3)

axes2[1].plot(semi_major_axis_scale, da_dt_pos, c='k')
axes2[1].plot(semi_major_axis_scale, da_dt_neg, c='k', alpha=0.3)

plt.show()