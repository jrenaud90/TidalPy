import TidalPy
import numpy as np
# Build the planet, let's choose Jupiter's moon Io
io = TidalPy.build_planet('Io')
io.paint()

jupiter = TidalPy.build_planet('Jupiter')
sol = TidalPy.build_planet('Sol')

# Change Io's orbital parameters if necessary
orbit = TidalPy.Orbit(sol, jupiter, io)

print("Io's orbital parameters:")
print(f"\te = {io.eccentricity}")
print(f"\tP = {orbit.rads2days(io.orbital_freq)} [days]")

# Let's make a helper function to quickly plot our findings
import matplotlib.pyplot as plt
def contour_plot(heating, rheology_name):
    fig, ax = plt.subplots()
    cpoints = [10, 11, 12, 13, 14, 15, 16]
    cinfo = ax.contourf(visco, shear, np.log10(heating), cpoints)
    ax.contour(visco, shear, np.log10(heating), cpoints, colors=('k',), linewidths=(0.5,))
    cbar = plt.colorbar(cinfo, ax=ax)
    cbar.ax.set_ylabel('Tidal Heating [W]')
    ax.set_ylabel('Shear Modulus [Pa]')
    ax.set_xlabel('Effective Viscosity [Pa s]')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_title(f'{rheology_name.title()} Rheology')
    plt.show()

# Load the tidal layer (mantle) with a pre-defined strength (Viscosity and Shear Modulus)
shear = np.logspace(4, 14, 100)
visco = np.logspace(8, 22, 100)
visco_mtx, shear_mtx = np.meshgrid(visco, shear)
io.mantle.set_strength(visco_mtx, shear_mtx)

# Pull out Tidal Heating and plot
tidal_heating = io.mantle.tidal_heating
contour_plot(tidal_heating, 'maxwell')

# Maxwell is the default rheology if none is provided in the planet's configuration
# But, we can quickly change the rheology without having to reload everything (or make a new config file)...
io.mantle.config['rheology'] = {'compliance_model': 'andrade'}
io.reinit()
tidal_heating = io.mantle.tidal_heating
contour_plot(tidal_heating, 'andrade')
