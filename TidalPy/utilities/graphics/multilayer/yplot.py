from typing import Union, List

import os
import numpy as np
import matplotlib.pyplot as plt

from ....exceptions import MissingArgumentError, IncorrectArgumentType

FILE_PATH = os.path.join(os.path.realpath(__file__), os.pardir)
PROP_CYCLE = plt.rcParams['axes.prop_cycle']
MP_COLORS = PROP_CYCLE.by_key()['color']

def yplot(tidal_ys: Union[List[np.ndarray], np.ndarray],
          radius : Union[List[np.ndarray], np.ndarray],
          depth_plot: bool = False, planet_radius: float = None,
          colors: List[str] = None, labels: List[str] = None,
          show_plot: bool = True, use_tobie_limits: bool = False,
          plot_tobie: bool = False, plot_roberts: bool = False):
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
    use_tobie_limits : bool = False
        If True, then plot will be made using Tobie2005 axis limits.

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

    fig_tidal_y, axes_tidal_y = plt.subplots(ncols=3, nrows=2, figsize=(7.5, 7.5))

    ax_y1 = axes_tidal_y[0, 0]
    ax_y2 = axes_tidal_y[0, 1]
    ax_y3 = axes_tidal_y[0, 2]
    ax_y4 = axes_tidal_y[1, 0]
    ax_y5 = axes_tidal_y[1, 1]
    ax_y6 = axes_tidal_y[1, 2]

    # TODO: Check units
    ax_y1.set(ylabel=ylabel, xlabel='$y_{1}$ [s$^{2}$ / m]', title='Radial Disp.')
    ax_y2.set(ylabel=ylabel, xlabel='$y_{2}$ [kg / m$^{3}$]', title='Radial Stress')
    ax_y3.set(ylabel=ylabel, xlabel='$y_{3}$ [s$^{2}$ / m]', title='Tang. Disp.')
    ax_y4.set(ylabel=ylabel, xlabel='$y_{4}$ [kg / m$^{3}$]', title='Tang. Stress.')
    ax_y5.set(ylabel=ylabel, xlabel='$y_{5}$ [unitless]', title='Grav. Potential Perturb.')
    ax_y6.set(ylabel=ylabel, xlabel='$y_{6}$ [1 / m]', title='Potential Stress')

    if not multiple_y:
        # Just throw them into a list to decrease amount of code.
        tidal_ys = [tidal_ys]
        radius = [radius]

    # Setup labels and colors
    if colors is None:
        colors = MP_COLORS

    if labels is None:
        labels = [f'y-{i}' for i in range(len(tidal_ys))]

    num_y = 0
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
        num_y += 1

    if multiple_y:
        if num_y < 4:
            ncols = 1
        elif num_y < 6:
            ncols = 2
        else:
            ncols = 3

        ax_y5.legend(ncol=ncols, fancybox=True, bbox_to_anchor=(1,-0.6), loc="lower right")

    if use_tobie_limits:
        ax_y1.set(xlim=(0.0, 0.15))
        ax_y2.set(xlim=(-2100, 4500))
        ax_y3.set(xlim=(-0.04, 0.04))
        ax_y4.set(xlim=(0, 2000))

    if plot_roberts:
        RN08_C = 'b'
        RN08_data = np.loadtxt(os.path.join(FILE_PATH, 'RN08-Data.csv'), skiprows=1, delimiter=',', dtype=str)
        RN08_data[RN08_data == ''] = np.nan
        RN08_data = RN08_data.astype(np.float)
        ax_y1.scatter(RN08_data[:, 0], RN08_data[:, 1] / 1000., label='RN08-HG', c=RN08_C, marker='.', s=20)
        ax_y1.scatter(RN08_data[:, 2], RN08_data[:, 3] / 1000., label='RN08-LC0', c=RN08_C, marker='+', s=20)

        ax_y3.scatter(RN08_data[:, 4], RN08_data[:, 5] / 1000., label='RN08-HG', c=RN08_C, marker='.', s=20)
        ax_y3.scatter(RN08_data[:, 6], RN08_data[:, 7] / 1000., label='RN08-LC0', c=RN08_C, marker='+', s=20)

        ax_y2.scatter(RN08_data[:, 8], RN08_data[:, 9] / 1000., label='RN08-HG', c=RN08_C, marker='.', s=20)
        ax_y2.scatter(RN08_data[:, 10], RN08_data[:, 11] / 1000., label='RN08-LC0', c=RN08_C, marker='+', s=20)

        ax_y4.scatter(RN08_data[:, 12], RN08_data[:, 13] / 1000., label='RN08-HG', c=RN08_C, marker='.', s=20)
        ax_y4.scatter(RN08_data[:, 14], RN08_data[:, 15] / 1000., label='RN08-LC0', c=RN08_C, marker='+', s=20)

    if plot_tobie:
        T05_C = 'r'
        T05_data = np.loadtxt(os.path.join(FILE_PATH, 'T05-Data.csv'), skiprows=1, delimiter=',', dtype=str)
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

    fig_tidal_y.tight_layout()
    if show_plot:
        plt.show()

    return fig_tidal_y
