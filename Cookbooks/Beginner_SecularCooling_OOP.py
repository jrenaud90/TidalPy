from TidalPy.planets import build_planet
from TidalPy.orbit import Orbit
from TidalPy.utilities.conversions import sec2myr, myr2sec
from TidalPy.utilities.progress import progress_bar

sol = build_planet('sol')
earth = build_planet('earth')
orbit = Orbit(sol, None, earth)

# Domain and initial conditions
initial_temperature = 1600.
time_domain_limits = (0., myr2sec(10000.))


# Equation to solve
@progress_bar(time_domain_limits)
def diffeq(time, variables, diffloop: bool = True):

    # Only one dependent variable in this example
    temperature = variables[0, :]

    # Load that temperature into the Earth's upper mantle, as well as time for any radiogenic calculations
    earth.time = sec2myr(time)
    earth.upper_mantle.temperature = temperature

    # Get out the temperature derivative and return it
    # TODO: Make this a state parameter so calculation occurs above?
    if diffloop:
        return earth.upper_mantle.calc_temperature_derivative()
    else:
        return earth.upper_mantle.radiogenic_heating

# Use Scipy package for integration
from scipy.integrate import solve_ivp
solution = solve_ivp(diffeq, time_domain_limits, [initial_temperature],
                     method='RK45', rtol=1e-6, vectorized=True, args=(True,))

solution_time = sec2myr(solution.t)
solution_temp = solution.y[0, :]

radiogenics = diffeq(solution.t, solution.y, diffloop=False)

# Plot the results
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
axtwin = ax.twinx()
ax.plot(solution_time, solution_temp, 'g')
axtwin.plot(solution_time, radiogenics/1e12, 'k')
ax.set_xlabel('Time [Myr]')
ax.set_ylabel('Upper Mantle Temperature [K]', c='g')
axtwin.set_ylabel('Radiogenic Heating [TW]')
axtwin.set_yscale('log')
plt.show()
