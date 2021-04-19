from typing import Union, List

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

from ....exceptions import MissingArgumentError, IncorrectArgumentType

def yplot(tidal_ys: Union[List[np.ndarray], np.ndarray],
          radius : Union[List[np.ndarray], np.ndarray],
          depth_plot: bool = False, planet_radius: float = None,
          colors: List[str] = None, labels: List[str] = None,
          show_plot: bool = True):
    """ Plot the six tidal y's in a six panel figure

    Parameters
    ----------
    tidal_ys : Union[List[np.ndarray], np.ndarray]
        Numpy array of tidal y's or list of numpy arrays of tidal y.
    radius : Union[List[np.ndarray], np.ndarray]
        Numpy array of radius or list of numpy arrays of radius.
    depth_plot : bool = False
        If True, will plot as a function of depth rather than radius.
    planet_radius : float = None
        Required if depth_plot is True.
    colors : List[str] = None
        Optional list of color names, one for each tidal_y solution.
    labels : List[str] = None
        Optional list of legend labels, one for each tidal_y solution.
    show_plot : bool = True
        If True, plt.show() will be called.

    Returns
    -------
    fig_tidal_y : plt.figure
    """

    ylabel = 'Radius [km]'
    if depth_plot:
        if planet_radius is None:
            MissingArgumentError('Planet radius is required for depth plot.')
        ylabel = 'Depth [km]'

    multiple_y = False
    if type(tidal_ys) == list:
        if type(radius) != list:
            IncorrectArgumentType('Radius must be provided as the same type as tidal_ys.')
        multiple_y = True

    fig_tidal_y, axes_tidal_y = plt.subplots(ncols=3, nrows=2, figsize=(10, 10))
    ax_y1 = axes_tidal_y[0, 0]
    ax_y2 = axes_tidal_y[0, 1]
    ax_y3 = axes_tidal_y[0, 2]
    ax_y4 = axes_tidal_y[1, 0]
    ax_y5 = axes_tidal_y[1, 1]
    ax_y6 = axes_tidal_y[1, 2]

    # TODO: Check units
    ax_y1.set(ylabel=ylabel, xlabel='$y_{1}$ [m / (m/s)$^{2}$]', title='Radial Disp.')
    ax_y2.set(ylabel=ylabel, xlabel='$y_{2}$ [kg / m$^{3}$', title='Radial Stress')
    ax_y3.set(ylabel=ylabel, xlabel='$y_{3}$ [m / (m/s)$^{2}$]', title='Tang. Disp.')
    ax_y4.set(ylabel=ylabel, xlabel='$y_{4}$ [kg / m$^{3}$]', title='Tang. Stress.')
    ax_y5.set(ylabel=ylabel, xlabel='$y_{5}$', title='Grav. Potential Perturb.')
    ax_y6.set(ylabel=ylabel, xlabel='$y_{6}$', title='Potential Stress')

    if not multiple_y:
        # Just throw them into a list to decrease amount of code.
        tidal_ys = [tidal_ys]
        radius = [radius]

    # Setup labels and colors
    if colors is None:
        colors = cm.get_cmap('jet')(np.linspace(0, 1, len(tidal_ys)))

    if labels is None:
        labels = [f'y-{i}' for i in range(len(tidal_ys))]

    for array_num, (radius_array, tidal_y_array) in enumerate(zip(radius, tidal_ys)):
        y1 = np.real(tidal_y_array[0, :])
        y2 = np.real(tidal_y_array[1, :])
        y3 = np.real(tidal_y_array[2, :])
        y4 = np.real(tidal_y_array[3, :])
        y5 = np.real(tidal_y_array[4, :])
        y6 = np.real(tidal_y_array[5, :])

        z = radius_array
        if depth_plot:
            z = planet_radius - radius_array

        ax_y1.plot(y1, z / 1000., label=labels[array_num], c=colors[array_num])
        ax_y2.plot(y2, z / 1000., label=labels[array_num], c=colors[array_num])
        ax_y3.plot(y3, z / 1000., label=labels[array_num], c=colors[array_num])
        ax_y4.plot(y4, z / 1000., label=labels[array_num], c=colors[array_num])
        ax_y5.plot(y5, z / 1000., label=labels[array_num], c=colors[array_num])
        ax_y6.plot(y6, z / 1000., label=labels[array_num], c=colors[array_num])

    if multiple_y:
        ax_y6.legend(loc='best')

    if show_plot:
        plt.show()

    return fig_tidal_y
