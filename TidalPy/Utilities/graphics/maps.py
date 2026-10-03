"""Draw surface maps of a field sampled on a colatitude-longitude grid, such as one radius of a world's
3D heating or stress-strain grid.

Cartopy supplies the projections when it is installed, matplotlib's Mollweide and rectangular axes
otherwise. No coastlines or other Natural Earth data are drawn, so nothing is downloaded.
"""

from __future__ import annotations

from typing import Dict, Optional, Tuple

import numpy as np
from matplotlib import colors as mpl_colors
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure

try:
    import cartopy.crs as ccrs
except ImportError:
    # Optional (the `graphics` extra); the matplotlib projections are used without it.
    ccrs = None

# Edit in place to restyle every later map.
MAP_PLOT_STYLE: Dict[str, object] = {
    "figure_size_inches": (9.0, 4.5),
    "colormap": "inferno",
    "symmetric_colormap": "RdBu_r",
    "grid_color": "0.6",
    "grid_alpha": 0.5,
    "grid_line_width": 0.5,
    "label_fontsize": 12,
    "title_fontsize": 14,
    "colorbar_fraction": 0.046,
    "colorbar_pad": 0.04,
}

# Projection name -> whether it needs cartopy.
MAP_PROJECTIONS: Dict[str, bool] = {
    "mollweide": False,
    "plate_carree": False,
    "robinson": True,
}


def resolve_map_backend(projection: str, central_longitude: float, use_cartopy: Optional[bool]) -> Tuple[str, bool]:
    """Validate a projection request and decide whether cartopy draws the map.

    Parameters
    ----------
    projection : str
        A key of `MAP_PROJECTIONS`, case-insensitive.
    central_longitude : float
        Longitude at the center of the map [deg].
    use_cartopy : bool or None
        True requires cartopy, False uses matplotlib, None uses cartopy when it is installed.

    Returns
    -------
    name : str
        The lower-case projection name.
    cartopy_used : bool
    """
    name = projection.lower()
    if name not in MAP_PROJECTIONS:
        raise ValueError(f"Unknown map projection '{projection}'; choose from {sorted(MAP_PROJECTIONS)}.")
    cartopy_installed = ccrs is not None
    cartopy_used = cartopy_installed if use_cartopy is None else bool(use_cartopy)
    if cartopy_used and not cartopy_installed:
        raise ImportError("Cartopy is not installed; install the `graphics` extra or pass use_cartopy=False.")
    if not cartopy_used:
        if MAP_PROJECTIONS[name]:
            raise ValueError(f"The '{name}' projection needs cartopy; use 'mollweide' or 'plate_carree' without it.")
        if central_longitude != 0.0:
            raise ValueError("Centering a map on a nonzero longitude needs cartopy.")
    return name, cartopy_used


def map_longitude_order(longitudes: np.ndarray, center: float) -> Tuple[np.ndarray, np.ndarray]:
    """Columns to draw and their longitudes, wrapped onto the 360 degrees around the map center.

    Wrapping around the center puts the seam of the data on the map edge. A longitude that repeats after
    wrapping, 0 and 2 pi say, keeps only its first column.
    
    Parameters
    ----------
    longitudes : array
        Longitude of each column [rad].
    center : float
        Longitude at the map center [deg].

    Returns
    -------
    column_order : numpy.ndarray of int
        Indices of the columns to draw, in increasing wrapped longitude.
    wrapped : numpy.ndarray of float
        Their wrapped longitudes [deg], within [center - 180, center + 180).
    """
    wrapped = (np.degrees(longitudes) - center + 180.0) % 360.0 - 180.0 + center
    column_order = np.argsort(wrapped, kind="stable")
    sorted_longitudes = wrapped[column_order]
    distinct = np.concatenate(([True], np.diff(sorted_longitudes) > 1.0e-9))
    return column_order[distinct], sorted_longitudes[distinct]


def map_cell_edges(centers: np.ndarray, lower: float, upper: float) -> np.ndarray:
    """Edges of the cells around sorted cell-centered samples, clipped to [lower, upper].

    Interior edges sit halfway between neighboring centers, the outer edges half a cell beyond the end
    centers. A single sample spans the whole range.
    """
    if centers.size == 1:
        return np.array([lower, upper])
    midpoints = 0.5 * (centers[1:] + centers[:-1])
    first = centers[0] - (midpoints[0] - centers[0])
    last = centers[-1] + (centers[-1] - midpoints[-1])
    return np.clip(np.concatenate(([first], midpoints, [last])), lower, upper)


def make_map_axes(
        nrows: int = 1,
        ncols: int = 1,
        projection: str = "mollweide",
        central_longitude: float = 0.0,
        use_cartopy: Optional[bool] = None,
        figure_size: Optional[Tuple[float, float]] = None,
        ) -> Tuple[Figure, np.ndarray]:
    """Create a figure of map panels for `plot_map`.

    Parameters
    ----------
    nrows, ncols : int, default 1
        Panel grid.
    projection : str, default "mollweide"
        A key of `MAP_PROJECTIONS`.
    central_longitude : float, default 0.0
        Longitude at the center of each panel [deg]; a nonzero value needs cartopy.
    use_cartopy : bool, optional
        True requires cartopy, False uses matplotlib; the default uses cartopy when it is installed.
    figure_size : (float, float), optional
        Figure size in inches; the default scales ``MAP_PLOT_STYLE["figure_size_inches"]`` by the panel grid.

    Returns
    -------
    figure : matplotlib.figure.Figure
    axes : numpy.ndarray of matplotlib.axes.Axes, shape (nrows, ncols)
        Cartopy GeoAxes when cartopy is used, else matplotlib Mollweide or rectangular axes.
    """
    name, cartopy_used = resolve_map_backend(projection, central_longitude, use_cartopy)
    if figure_size is None:
        width, height = MAP_PLOT_STYLE["figure_size_inches"]
        figure_size = (width * ncols, height * nrows)
    figure = plt.figure(figsize=figure_size)
    axes = np.empty((nrows, ncols), dtype=object)
    for index in range(nrows * ncols):
        if cartopy_used:
            if name == "mollweide":
                panel_projection = ccrs.Mollweide(central_longitude=central_longitude)
            elif name == "robinson":
                panel_projection = ccrs.Robinson(central_longitude=central_longitude)
            else:
                panel_projection = ccrs.PlateCarree(central_longitude=central_longitude)
            axis = figure.add_subplot(nrows, ncols, index + 1, projection=panel_projection)
        elif name == "mollweide":
            axis = figure.add_subplot(nrows, ncols, index + 1, projection="mollweide")
        else:
            axis = figure.add_subplot(nrows, ncols, index + 1)
        axes.flat[index] = axis
    return figure, axes


def plot_map(
        longitudes: np.ndarray,
        colatitudes: np.ndarray,
        values: np.ndarray,
        title: Optional[str] = None,
        colorbar_label: Optional[str] = None,
        colormap: Optional[str] = None,
        projection: str = "mollweide",
        central_longitude: float = 0.0,
        value_limits: Optional[Tuple[float, float]] = None,
        symmetric: bool = False,
        log_scale: bool = False,
        grid_lines: bool = True,
        colorbar: bool = True,
        axis: Optional[Axes] = None,
        use_cartopy: Optional[bool] = None,
        show_plot: bool = False,
        ) -> Tuple[Figure, Axes]:
    """Draw a field sampled on a colatitude-longitude grid as a global map.

    Parameters
    ----------
    longitudes : array
        Longitude of each column [rad], wrapped onto the 360 degrees around the map center and sorted.
    colatitudes : array
        Colatitude of each row [rad], 0 at the north pole.
    values : array, shape (len(colatitudes), len(longitudes))
        The field, in the axis order of the world 3D grids. NaN cells are left blank.
    title, colorbar_label : str, optional
        Panel title and colorbar label.
    colormap : str, optional
        Matplotlib colormap; defaults to ``MAP_PLOT_STYLE["symmetric_colormap"]`` when `symmetric` and to
        ``MAP_PLOT_STYLE["colormap"]`` otherwise.
    projection : str, default "mollweide"
        A key of `MAP_PROJECTIONS`; ignored when `axis` is given.
    central_longitude : float, default 0.0
        Longitude at the map center [deg]; a nonzero value needs cartopy. Ignored when `axis` is given,
        which already has a center.
    value_limits : (float, float), optional
        Color scale limits; the default spans the finite values.
    symmetric : bool, default False
        Center the color scale on zero, for signed fields such as stress and strain.
    log_scale : bool, default False
        Logarithmic color scale over the positive values, for fields such as heating; other cells are left blank.
    grid_lines : bool, default True
        Draw latitude and longitude lines.
    colorbar : bool, default True
        Add a colorbar beside the panel.
    axis : matplotlib.axes.Axes, optional
        Panel to draw into, from `make_map_axes`; a new single-panel figure is made when omitted.
    use_cartopy : bool, optional
        True requires cartopy, False uses matplotlib; the default uses cartopy when it is installed. Ignored when
        `axis` is given.
    show_plot : bool, default False
        Call ``matplotlib.pyplot.show()`` before returning.

    Returns
    -------
    figure : matplotlib.figure.Figure
    axis : matplotlib.axes.Axes

    Raises
    ------
    ValueError
        `values` does not match the axes, has nothing to draw, the color scale options conflict, or the
        backend cannot draw the requested projection.
    ImportError
        `use_cartopy` is True and cartopy is not installed.
    """
    style = MAP_PLOT_STYLE
    longitude_values = np.asarray(longitudes, dtype=np.float64).ravel()
    colatitude_values = np.asarray(colatitudes, dtype=np.float64).ravel()
    field = np.asarray(values, dtype=np.float64)
    if field.shape != (colatitude_values.size, longitude_values.size):
        raise ValueError(
            f"`values` must be shaped (len(colatitudes), len(longitudes)) = "
            f"({colatitude_values.size}, {longitude_values.size}); got {field.shape}.")
    if not (np.all(np.isfinite(longitude_values)) and np.all(np.isfinite(colatitude_values))):
        raise ValueError("`longitudes` and `colatitudes` must be finite.")
    if symmetric and log_scale:
        raise ValueError("`symmetric` and `log_scale` cannot be combined.")
    drawable = np.isfinite(field) & ((field > 0.0) if log_scale else True)
    if not np.any(drawable):
        raise ValueError("`values` has no finite value to draw" + (" above zero on a log scale." if log_scale else "."))

    if axis is None:
        figure, axes = make_map_axes(
            projection=projection,
            central_longitude=central_longitude,
            use_cartopy=use_cartopy)
        axis = axes[0, 0]
    else:
        figure = axis.figure
    cartopy_axis = ccrs is not None and isinstance(getattr(axis, "projection", None), ccrs.Projection)
    center = float(axis.projection.proj4_params.get("lon_0", 0.0)) if cartopy_axis else 0.0

    # Columns in increasing longitude around the map center, rows in increasing latitude.
    column_order, longitude_degrees = map_longitude_order(longitude_values, center)
    latitude_degrees = 90.0 - np.degrees(colatitude_values)
    row_order = np.argsort(latitude_degrees, kind="stable")
    latitude_degrees = latitude_degrees[row_order]
    field = field[np.ix_(row_order, column_order)]

    masked = np.ma.masked_invalid(field)
    if log_scale:
        masked = np.ma.masked_less_equal(masked, 0.0)
    drawn = masked.compressed()
    if value_limits is not None:
        lower, upper = value_limits
    elif symmetric:
        upper = float(np.max(np.abs(drawn)))
        lower = -upper
    else:
        lower, upper = float(drawn.min()), float(drawn.max())
    norm = mpl_colors.LogNorm(vmin=lower, vmax=upper) if log_scale else mpl_colors.Normalize(vmin=lower, vmax=upper)
    if colormap is None:
        colormap = style["symmetric_colormap"] if symmetric else style["colormap"]

    longitude_edges = map_cell_edges(longitude_degrees, center - 180.0, center + 180.0)
    latitude_edges = map_cell_edges(latitude_degrees, -90.0, 90.0)
    grid_style = dict(color=style["grid_color"], alpha=style["grid_alpha"], linewidth=style["grid_line_width"])
    if cartopy_axis:
        mesh = axis.pcolormesh(
            longitude_edges,
            latitude_edges,
            masked,
            cmap=colormap,
            norm=norm,
            shading="flat",
            transform=ccrs.PlateCarree())
        axis.set_global()
        if grid_lines:
            axis.gridlines(**grid_style)
    elif axis.name == "mollweide":
        # The matplotlib geographic axes take radians and label longitudes across the equator, so fewer.
        mesh = axis.pcolormesh(
            np.radians(longitude_edges),
            np.radians(latitude_edges),
            masked,
            cmap=colormap,
            norm=norm,
            shading="flat")
        axis.set_xticks(np.radians(np.arange(-120.0, 121.0, 60.0)))
        axis.grid(grid_lines, **grid_style)
    else:
        mesh = axis.pcolormesh(
            longitude_edges,
            latitude_edges,
            masked,
            cmap=colormap,
            norm=norm,
            shading="flat")
        axis.set_xlim(-180.0, 180.0)
        axis.set_ylim(-90.0, 90.0)
        axis.set_xlabel("Longitude [deg]", fontsize=style["label_fontsize"])
        axis.set_ylabel("Latitude [deg]", fontsize=style["label_fontsize"])
        axis.grid(grid_lines, **grid_style)

    if colorbar:
        figure.colorbar(
            mesh,
            ax=axis,
            fraction=style["colorbar_fraction"],
            pad=style["colorbar_pad"],
            label=colorbar_label)
    if title is not None:
        axis.set_title(title, fontsize=style["title_fontsize"])
    if show_plot:
        plt.show()
    return figure, axis
