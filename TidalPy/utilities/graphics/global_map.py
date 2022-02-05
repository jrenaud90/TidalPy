""" Functions to project data onto a 2D representation of a sphere. """
import os.path
import numpy as np

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from typing import Union, List
from .helper import get_cmap
from ...io_helper import unique_path

ROTATED_POLE_DEFAULT_INPUT = {
    'pole_latitude': 45,
    'pole_longitude': 180
    }

_KNOWN_PROJECTIONS = {
    'PlateCarree': ccrs.PlateCarree,
    'AlbersEqualArea': ccrs.AlbersEqualArea,
    'AzimuthalEquidistant': ccrs.AzimuthalEquidistant,
    'EquidistantConic': ccrs.EquidistantConic,
    'LambertConformal': ccrs.LambertConformal,
    'LambertCylindrical': ccrs.LambertCylindrical,
    'Mercator': ccrs.Mercator,
    'Miller': ccrs.Miller,
    'Mollweide': ccrs.Mollweide,
        # A Mollweide projection. This projection is pseudocylindrical, and equal area.
        # Parallels are unequally-spaced straight lines, while meridians are elliptical arcs up to
        # semicircles on the edges. Poles are points.
        # It is commonly used for world maps, or interrupted with several central meridians.
    'Orthographic': ccrs.Orthographic,
    'Robinson': ccrs.Robinson,
        # A Robinson projection. This projection is pseudocylindrical, and a compromise
        # that is neither equal-area nor conformal. Parallels are unequally-spaced straight lines,
        # and meridians are curved lines of no particular form. It is commonly used for “visually-appealing” world maps.
    'Sinusoidal': ccrs.Sinusoidal,
        # A Sinusoidal projection. This projection is equal-area.
    'RotatedPole': ccrs.RotatedPole,
        # A rotated latitude/longitude projected coordinate system with cylindrical topology and projected distance.
        # Coordinates are measured in projection metres.
    'NearsidePerspective': ccrs.NearsidePerspective,
        # Perspective view looking directly down from above a point on the globe.
        # In this projection, the projected coordinates are x and y measured from the origin of a plane tangent to
        # the Earth directly below the perspective point (e.g. a satellite).
    }

# Make a copy of the projection keys as lower case to avoid key errors.
KNOWN_PROJECTIONS = dict()
for _projection_name, _projection in _KNOWN_PROJECTIONS.items():
    KNOWN_PROJECTIONS[_projection_name] = _projection
    if _projection_name.lower() not in KNOWN_PROJECTIONS:
        KNOWN_PROJECTIONS[_projection_name.lower()] = _projection

def projection_map(longitude: np.ndarray, colatitude: np.ndarray, data: np.ndarray,
                   cpoints: Union[int, List[float]] = 30,
                   figure_scale: float = 1., aspect_ratio: float = 2.,
                   cmap: str = 'vik',
                   xlabel: str = 'Longitude [deg.]', ylabel: str = 'Latitude [deg.]',
                   zlabel: str = '', title : str = '', zlog: bool = False,
                   zticks: list = None, ztick_labels: list = None,
                   projection: str = 'Mollweide', original_data_projection: str = 'PlateCarree',
                   show_earth_coast: bool = False, show_grid_lines: bool = True,
                   auto_save: bool = True, auto_show: bool = True, save_png: bool = False,
                   png_dpi: int = 300, filename: str = 'unknown_projection_plot', rotated_pole_input: dict = None):
    """

    Parameters
    ----------
    longitude
    colatitude
    data
    projection
    add_earth_coast
    add_latlong_lines
    figure_scale
    aspect_ratio

    Returns
    -------

        # @manual{Cartopy,
        # author = {{Met Office}},
        # title = {Cartopy: a cartographic python library with a Matplotlib interface},
        # year = {2010 - 2015},
        # address = {Exeter, Devon },
        # url = {https://scitools.org.uk/cartopy},
        # doi = {10.5281/zenodo.1182735}
        # }
    """

    # Find projection
    if projection.lower() not in KNOWN_PROJECTIONS:
        raise KeyError(f'Unknown projection model for project_map: {projection}')
    projection = KNOWN_PROJECTIONS[projection.lower()]

    # Find data projection
    if original_data_projection.lower() not in KNOWN_PROJECTIONS:
        raise KeyError(f'Unknown data projection model for project_map: {projection}')
    data_projection = KNOWN_PROJECTIONS[original_data_projection.lower()]

    # Check rotated pole data
    if projection is ccrs.RotatedPole:
        if rotated_pole_input is None:
            rotated_pole_input = ROTATED_POLE_DEFAULT_INPUT
        else:
            rotated_pole_input = {**ROTATED_POLE_DEFAULT_INPUT, **rotated_pole_input}
        projection_instance = projection(**rotated_pole_input)
    else:
        rotated_pole_input = None
        projection_instance = projection()

    # Find color map
    cmap = get_cmap(cmap)

    # Make figure
    fig = plt.figure(figsize=(figure_scale * aspect_ratio * 4., figure_scale * 4.))
    ax = plt.axes(projection=projection_instance)

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
    cbdata = ax.contourf(longitude, latitude, data.T, cpoints, cmap=cmap,
                         transform=data_projection())

    # Make colorbar
    colorbar = plt.colorbar(cbdata, ax=ax)
    if zticks is not None:
        colorbar.set_ticks(zticks)
    if ztick_labels is not None:
        colorbar.ax.set_yticklabels(ztick_labels)

    # Plot lat and long gridlines
    if show_grid_lines:
        ax.gridlines()

    # Set labels
    ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
    colorbar.set_label(zlabel)

    # Save figure
    if auto_save:
        if '.pdf' == filename[-4:]:
            filename = filename[:-4]

        filename = unique_path(filename + '.pdf', is_dir=False)
        filepath = filename[:-4]
        fig.savefig(filepath + '.pdf')
        if save_png:
            fig.savefig(filepath + '.png', dpi=png_dpi)

    # Show figure
    if auto_show:
        plt.show()

    return fig, ax
