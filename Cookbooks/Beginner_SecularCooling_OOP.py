from TidalPy.planets import build_planet
from TidalPy.orbit import Orbit
from TidalPy.utilities.conversions import sec2myr, myr2sec
from TidalPy.utilities.progress import progress_bar

sol = build_planet('sol')
earth = build_planet('earth')
orbit = Orbit(sol, None, earth)

# Domain and initial conditions
initial_temperature_upper = 1600.
initial_temperature_lower = 1900.
initial_conditions = [initial_temperature_lower, initial_temperature_upper]

# Uncomment the below to see the Earth's upper mantle adjust to the lower mantle and initial conditions.
time_domain_limits = (0., myr2sec(50))

# Uncomment the below to see the Earth cool to the surface temperature.
# time_domain_limits = (0., myr2sec(5e5))

# Uncomment the below lines to see the affect of a non-zero molar activation volume (pressure-dependent viscosity).
# time_domain_limits = (0., myr2sec(2500))
# earth.config['layers']['Lower_Mantle']['rheology']['solid_viscosity']['molar_activation_volume'] = -1.5e-6
# earth.reinit()


# Equation to solve
@progress_bar(time_domain_limits)
def diffeq(time, variables, diff_loop: bool = True):

    # Only one dependent variable in this example
    lower_temperature = variables[0, :]
    upper_temperature = variables[1, :]

    # Load that temperature into the Earth's upper mantle, as well as time for any radiogenic calculations
    earth.time = sec2myr(time)

    # The current OOP methodology propagates temperature *down* (the lower layers need to know their upper boundary
    #  condition in order to calculate cooling). Therefore, start with the uppermost layers.
    earth.upper_mantle.temperature = upper_temperature
    earth.lower_mantle.temperature = lower_temperature

    # Get out the temperature derivative and return it
    # TODO: Make this a state parameter so calculation occurs above?
    if diff_loop:
        lower_dTdt = earth.lower_mantle.calc_temperature_derivative()
        upper_dTdt = earth.upper_mantle.calc_temperature_derivative()
        return [lower_dTdt, upper_dTdt]
    else:
        lower_radio = earth.lower_mantle.radiogenic_heating
        upper_radio = earth.upper_mantle.radiogenic_heating
        return upper_radio, lower_radio

# Use Scipy package for integration
from scipy.integrate import solve_ivp
solution = solve_ivp(diffeq, time_domain_limits, initial_conditions,
                     method='LSODA', rtol=1e-6, vectorized=True)
print() # adds a space after the progress bar

solution_time = sec2myr(solution.t)
solution_temp_lower = solution.y[0, :]
solution_temp_upper = solution.y[1, :]

radiogenics_lower, radiogenics_upper = diffeq(solution.t, solution.y, diff_loop=False)

# Plot the results
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
axtwin = ax.twinx()
axtwin.plot(solution_time, radiogenics_lower/1e12, 'k-', label='Lower Mantle')
axtwin.plot(solution_time, radiogenics_upper/1e12, 'k--', label='Upper Mantle')
ax.plot(solution_time, solution_temp_lower, 'g-', label='Lower Mantle')
ax.plot(solution_time, solution_temp_upper, 'g--', label='Upper Mantle')
ax.set_xlabel('Time [Myr]')
ax.set_ylabel('Temperature [K]', c='g')
axtwin.set_ylabel('Radiogenic Heating [TW]')
ax.legend(loc='best')
plt.show()
