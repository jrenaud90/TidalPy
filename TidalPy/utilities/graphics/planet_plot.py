import numpy as np
from matplotlib import gridspec as gspec, pyplot as plt

from ...exceptions import ParameterMissingError

SCATTER_SIZE = 1
SCALE = 3


def geotherm_plot(radii: np.ndarray,
                  gravitys: np.ndarray, pressures: np.ndarray, densitys: np.ndarray, temperatures: np.ndarray = None,
                  planet_radius: float = None, bulk_density: float = None, planet_name: str = None,
                  depth_plot: bool = False, auto_show: bool = False, annotate: bool = True):
    """ Plots the depth plot of a planet in 3 or 4 panels (temperature is optional)

    Parameters
    ----------
    radii : np.ndarray
        Planet's radius discretized into multiple chunks in a 1-D array [m]
    gravitys : np.ndarray
        Planet's gravity as a function of depth or radius
    pressures : np.ndarray
        Planet's pressure as a function of depth or radius
    densitys : np.ndarray
        Planet's density as a function of depth or radius
    temperatures : np.ndarray
        (optional) Planet's internal temperature as a function of depth or radius
    planet_radius : float
        (optional) Planet's outer radius - used for depth plot
    bulk_density : float
        (optional) Planet's bulk density can be shown as an annotation
    planet_name : str
        Planet's name
    depth_plot : bool
        (default = False) If True then x-axis will be depth not radius
    auto_show : bool
        (default = False) If True then plt.show will be called
    annotate : bool
        (default = True) Adds annotations to some of the plots (surface gravity, base pressure, etc.)

    Returns
    -------
    fig : matplotlib.Figure

    """

    # Determine if this is a depth or radius plot
    if depth_plot:
        if planet_radius is None:
            raise ParameterMissingError
        x = planet_radius - radii
    else:
        x = radii
    shape = x.shape

    assert gravitys.shape == shape
    assert pressures.shape == shape
    assert densitys.shape == shape
    use_temperature = False
    if temperatures is not None:
        assert temperatures.shape == shape
        use_temperature = True

    # Perform any unit conversions
    pressures = pressures / 1e9  # Use GPa
    densitys = densitys / 1000
    x = x / 1000  # Use km
    if bulk_density is not None:
        bulk_density = bulk_density / 1000

    # If depth_plot then the plots will be stacked vertically, otherwise they will be arranged horizontally.
    n = 3
    if use_temperature:
        n = 4

    # Construct plot grid.
    fig = plt.figure(figsize=(n * SCALE, 1 * SCALE))
    gs = gspec.GridSpec(nrows=1, ncols=n, wspace=0.4, hspace=0.2, figure=fig)
    axes = [fig.add_subplot(gs[i]) for i in range(n)]
    ys = [gravitys, pressures, densitys]
    colors = ['g', 'b', 'k']
    names = ['Gravity [m s$^{2}$]', 'Pressure [GPa]', 'Density [10$^{3}$ kg m$^{-3}$]']
    if use_temperature:
        ys.append(temperatures)
        colors.append('r')
        names.append('Temperature [K]')

    for plot_i, (y, ax, color, name) in enumerate(zip(ys, axes, colors, names)):

        x_ = x
        y_ = y
        if depth_plot:
            # Depth plot has depth in the y axis and it is reversed
            x_ = y
            y_ = x
            ax.set_ylim(max(y_), 0)

        # Main plot
        ax.scatter(x_, y_, s=SCATTER_SIZE, c=color)

        # Plot style
        if depth_plot:
            ax.set_ylabel('Depth [km]')
            ax.set_xlabel(name)
        else:
            ax.set_ylabel(name)
            ax.set_xlabel('Radius [km]')

        # Add annotations
        if annotate:
            if y is gravitys:
                if depth_plot:
                    pos = (0.05, 0.9)
                else:
                    pos = (0.05, 0.9)
                ax.text(*pos, f'surf = {gravitys[-1]:0.2f}',
                        horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
            if y is pressures:
                if depth_plot:
                    pos = (0.05, 0.15)
                else:
                    pos = (0.05, 0.15)
                ax.text(*pos, f'base = {pressures[0]:0.2f}',
                        horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
            if y is densitys and bulk_density is not None:
                if depth_plot:
                    pos = (0.05, 0.15)
                else:
                    pos = (0.05, 0.15)
                ax.text(*pos, f'bulk = {bulk_density:0.2f}',
                        horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)

    if planet_name is not None:
        axes[0].set_title(planet_name)

    gs.tight_layout(fig)
    if auto_show:
        plt.show()

    return fig
