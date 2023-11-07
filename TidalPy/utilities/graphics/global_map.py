""" Functions to project data onto a 2D representation of a sphere. """
from typing import List, Union
import warnings

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Colormap

from TidalPy.paths import unique_path

cartopy_installed = False
try:
    import cartopy.crs as ccrs
    cartopy_installed = True
    _KNOWN_PROJECTIONS = {
        # # Matplotlib Projections
        'rectilinear'         : None,

        # # CCR Projections
        'PlateCarree'         : ccrs.PlateCarree,
        'AlbersEqualArea'     : ccrs.AlbersEqualArea,
        'AzimuthalEquidistant': ccrs.AzimuthalEquidistant,
        'EquidistantConic'    : ccrs.EquidistantConic,
        'LambertConformal'    : ccrs.LambertConformal,
        'LambertCylindrical'  : ccrs.LambertCylindrical,
        'Mercator'            : ccrs.Mercator,
        'Miller'              : ccrs.Miller,
        'Mollweide'           : ccrs.Mollweide,
        # A Mollweide projection. This projection is pseudocylindrical, and equal area.
        # Parallels are unequally-spaced straight lines, while meridians are elliptical arcs up to
        # semicircles on the edges. Poles are points.
        # It is commonly used for world maps, or interrupted with several central meridians.
        'Orthographic'        : ccrs.Orthographic,
        'Robinson'            : ccrs.Robinson,
        # A Robinson projection. This projection is pseudocylindrical, and a compromise
        # that is neither equal-area nor conformal. Parallels are unequally-spaced straight lines,
        # and meridians are curved lines of no particular form. It is commonly used for “visually-appealing” world maps.
        'Sinusoidal'          : ccrs.Sinusoidal,
        # A Sinusoidal projection. This projection is equal-area.
        'RotatedPole'         : ccrs.RotatedPole,
        # A rotated latitude/longitude projected coordinate system with cylindrical topology and projected distance.
        # Coordinates are measured in projection metres.
        'NearsidePerspective' : ccrs.NearsidePerspective,
        # Perspective view looking directly down from above a point on the globe.
        # In this projection, the projected coordinates are x and y measured from the origin of a plane tangent to
        # the Earth directly below the perspective point (e.g. a satellite).
        }
except ImportError:
    # Cartopy install not found.
    warnings.warn('The Cartopy package can not be found. 2D projection maps will not be available. Install cartopy '
                  'by using the Anaconda distribution of Python and then calling `conda install cartopy`.')
    class ccrs:
        pass

    _KNOWN_PROJECTIONS = {
        # # Matplotlib Projections
        'rectilinear': None
        }

ROTATED_POLE_DEFAULT_INPUT = {
    'pole_latitude' : 45,
    'pole_longitude': 180
    }

# Make a copy of the projection keys as lower case to avoid key errors.
KNOWN_PROJECTIONS = dict()
for _projection_name, _projection in _KNOWN_PROJECTIONS.items():
    KNOWN_PROJECTIONS[_projection_name] = _projection
    if _projection_name.lower() not in KNOWN_PROJECTIONS:
        KNOWN_PROJECTIONS[_projection_name.lower()] = _projection


def projection_map(
    longitude: np.ndarray, colatitude: np.ndarray, data: np.ndarray,
    cpoints: Union[int, List[float]] = 30,
    figure_scale: float = 1., aspect_ratio: float = 2.,
    cmap: Union[str, Colormap] = 'Blues',
    xlabel: str = 'Longitude [deg.]', ylabel: str = 'Latitude [deg.]',
    zlabel: str = '', title: str = '', zlog: bool = False,
    zticks: list = None, ztick_labels: list = None, horizontal_cbar: bool = False,
    projection: str = 'Mollweide', original_data_projection: str = 'PlateCarree',
    show_earth_coast: bool = False, show_grid_lines: bool = True,
    ax: plt.Axes = None, cbax: plt.Axes = None,
    auto_save: bool = False, auto_show: bool = True, save_png: bool = False,
    png_dpi: int = 300, filename: str = 'unknown_projection_plot', rotated_pole_input: dict = None,
    central_longitude: float = 0., dark_mode: bool = False
    ):
    """ Quickly plot surface maps using various kinds of latitude and longitude projections.

    This tool is made possible by the awesome work of the `Cartopy` package. If you use this feature please consider
    citing this package:
    [Cartopy Website](https://scitools.org.uk/cartopy/docs/latest/index.html)
    DOI: [10.5281/zenodo.1182735](https://doi.org/10.5281/zenodo.5842769)

    Projection Options
    ------------------
    'PlateCarree': ccrs.PlateCarree,
    'AlbersEqualArea': ccrs.AlbersEqualArea,
    'AzimuthalEquidistant': ccrs.AzimuthalEquidistant,
    'EquidistantConic': ccrs.EquidistantConic,
    'LambertConformal': ccrs.LambertConformal,
    'LambertCylindrical': ccrs.LambertCylindrical,
    'Mercator': ccrs.Mercator,
    'Miller': ccrs.Miller,
    'Mollweide': ccrs.Mollweide,
        A Mollweide projection. This projection is pseudocylindrical, and equal area.
        Parallels are unequally-spaced straight lines, while meridians are elliptical arcs up to
        semicircles on the edges. Poles are points.
        It is commonly used for world maps, or interrupted with several central meridians.
    'Orthographic': ccrs.Orthographic,
    'Robinson': ccrs.Robinson,
        A Robinson projection. This projection is pseudocylindrical, and a compromise
        that is neither equal-area nor conformal. Parallels are unequally-spaced straight lines,
        and meridians are curved lines of no particular form. It is commonly used for “visually-appealing” world maps.
    'Sinusoidal': ccrs.Sinusoidal,
        A Sinusoidal projection. This projection is equal-area.
    'RotatedPole': ccrs.RotatedPole,
        A rotated latitude/longitude projected coordinate system with cylindrical topology and projected distance.
        Coordinates are measured in projection metres.
    'NearsidePerspective': ccrs.NearsidePerspective,
        Perspective view looking directly down from above a point on the globe.
        In this projection, the projected coordinates are x and y measured from the origin of a plane tangent to
        the Earth directly below the perspective point (e.g. a satellite).

    Parameters
    ----------
    longitude : np.ndarray
        Longitude at which `data` was calculated [degs].
    colatitude : np.ndarray
        Colatitude at which `data` was calculated [degs].
    data : np.ndarray
        2D Data calculated at each longitude and latitude.
    cpoints : Union[int, List[float]] = 30
        Contour plot points.
    figure_scale : int = 1.
        Scale of the matplotlib figure.
    aspect_ratio : float = 2.
        Aspect ratio (w/h) of the matplotlib figure.
    cmap : Union[str, Colormap] = 'vik'
        Colormap used for the contour plot.
    xlabel : str = ''
        X-axis label.
    ylabel : str = ''
        Y-axis label.
    zlabel : str = ''
        Color bar label.
    title : str = ''
        Figure title.
    zlog : bool = False
        If True, `data` will be log-scaled.
    zticks : list = None
        If provided, these ticks will be used for the colorbar over matplotlib's default.
    ztick_labels : list = None
        If provided, these tick labels will be used for the colorbar over matplotlib's default.
    horizontal_cbar: bool = False
        If True, then the colorbar will be placed underneath the plot, horizontally.
    projection : str = 'Mollweide'
        Desired map projection. See function description for options.
    original_data_projection : str = 'PlateCarree'
        The projection that the original `data` was calculated using. For evenly spaced grid use 'PlateCarree'.
    show_earth_coast : bool = False
        If True, the modern-day Earth's coastline will be plotted over top of the data. This can be useful for
        comparison purposes.
    show_grid_lines : bool = True
        If True, latitude and longitude lines will be plotted over the data.
    ax : plt.Axes = None
        If a matplotlib figure axis is provided then the function will use it rather than creating a new figure.
        Note, this disables saving and show figure.
    cbax : plt.Axes = None
        Colorbar axis
    auto_save : bool = True
        If True, the figure will be saved to disk.
    auto_show : bool = True
        If True, the figure will be shown using `plt.show()`.
    save_png : bool = False
        If True, then a .png image will be saved alongside the standard pdf.
    png_dpi : int = 300
        PNG dots per inch.
    filename : str = 'unknown_projection_plot'
        Save name for the image. Can include a directory path.
    rotated_pole_input : dict = None
        Optional inputs for the Rotated Pole projection.
    central_longitude: float = 0.
        Longitude at the center of the plot.

    Returns
    -------
    fig: plt.Figure
        Image figure.
    ax: plt.Axis
        Image matplotlib axis.

    """

    # Dark mode
    if dark_mode:
        plt.style.use('dark_background')

    # Find projection
    if projection.lower() not in KNOWN_PROJECTIONS:
        raise KeyError(f'Unknown projection model for project_map: {projection}')
    projection = KNOWN_PROJECTIONS[projection.lower()]

    # Find data projection
    if original_data_projection.lower() not in KNOWN_PROJECTIONS:
        raise KeyError(f'Unknown data projection model for project_map: {projection}')
    data_projection = KNOWN_PROJECTIONS[original_data_projection.lower()]

    # Check rotated pole data
    projection_instance = None
    if cartopy_installed:
        if projection is ccrs.RotatedPole:
            if rotated_pole_input is None:
                rotated_pole_input = ROTATED_POLE_DEFAULT_INPUT
            else:
                rotated_pole_input = {**ROTATED_POLE_DEFAULT_INPUT, **rotated_pole_input}
            try:
                # Some versions of cartopy do not have the central longitude argument.
                projection_instance = projection(central_longitude=central_longitude, **rotated_pole_input)
            except TypeError:
                projection_instance = projection(**rotated_pole_input)
    if projection_instance is None:
        if projection is None:
            rotated_pole_input = None
            projection_instance = None
        else:
            rotated_pole_input = None
            projection_instance = projection(central_longitude=central_longitude)

    # Make figure
    if ax is None:
        fig = plt.figure(figsize=(figure_scale * aspect_ratio * 4., figure_scale * 4.))
        if projection_instance is None:
            ax = plt.axes(projection=None)
        else:
            ax = plt.axes(projection=projection_instance)
        premade_ax = False
    else:
        fig = None
        premade_ax = True

    # make the map global rather than have it zoom in to the extents of any plotted data
    ax.set_global()

    # Add in an Earth-like coastline for comparison purposes.
    if show_earth_coast:
        ax.coastlines()

    # We need to plot in latitude instead of colatitude.
    latitude = 90 - colatitude

    # Log data if needed.
    if zlog:
        data = np.log10(data)

    # Plot data making sure to transform it to the correct coordinate system.
    if data_projection is None:
        cbdata = ax.contourf(
            longitude, latitude, data.T, cpoints, cmap=cmap,
            transform=None
            )
    else:
        cbdata = ax.contourf(
            longitude, latitude, data.T, cpoints, cmap=cmap,
            transform=data_projection()
            )

    # Make colorbar
    if horizontal_cbar:
        location = 'bottom'
    else:
        location = 'right'

    if cbax is None:
        colorbar = plt.colorbar(cbdata, ax=ax, location=location)
    else:
        colorbar = plt.colorbar(cbdata, cax=cbax)
    if zticks is not None:
        colorbar.set_ticks(zticks)
    if ztick_labels is not None:
        colorbar.ax.set_yticklabels(ztick_labels)

    # Set labels
    ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
    colorbar.set_label(zlabel)

    # Plot lat and long gridlines
    if show_grid_lines:
        gl = ax.gridlines(
            draw_labels={'left':True, 'bottom':False, 'geo':True}, linestyle='-', alpha=0.35, x_inline=False,
            xlocs=[-120, -60, 0, 60, 120],
            ylocs=[-60, -30, 0, 30, 60]
            )
        # gl.rotate_labels = False
        if dark_mode:
            gl.xlabel_style = {'color': 'black'}
        else:
            gl.xlabel_style = {'color': 'white'}
        gl.right_labels = False
        gl.top_labels = False
        gl.bottom_labels = False


    # Save figure
    if auto_save and not premade_ax:
        if '.pdf' == filename[-4:]:
            filename = filename[:-4]

        filename = unique_path(filename + '.pdf', is_dir=False)
        filepath = filename[:-4]
        fig.savefig(filepath + '.pdf')
        if save_png:
            fig.savefig(filepath + '.png', dpi=png_dpi)

    # Show figure
    if auto_show and not premade_ax:
        plt.show()

    return fig, ax
