"""Tests for the surface map helpers ``plot_map`` and ``make_map_axes`` in TidalPy.Utilities.graphics."""
import matplotlib

matplotlib.use("Agg")

import numpy as np
import pytest
from matplotlib import pyplot as plt

from TidalPy.Utilities.graphics import MAP_PLOT_STYLE, MAP_PROJECTIONS, make_map_axes, plot_map
from TidalPy.Utilities.graphics import maps as maps_module

# A longitude grid that repeats 0 and 2 pi, as np.linspace(0, 2 pi, n) does, and a colatitude grid off the poles.
LONGITUDES = np.linspace(0.0, 2.0 * np.pi, 13)
COLATITUDES = np.linspace(0.1, np.pi - 0.1, 7)
COLATITUDE_GRID, LONGITUDE_GRID = np.meshgrid(COLATITUDES, LONGITUDES, indexing="ij")
FIELD = np.sin(COLATITUDE_GRID) ** 2 * np.cos(2.0 * LONGITUDE_GRID) + 0.1 * np.cos(COLATITUDE_GRID)


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close("all")


def expected_drawn_field(longitudes, colatitudes, values):
    """The field in drawing order: rows by increasing latitude, columns by longitude wrapped onto [-180, 180)."""
    wrapped = (np.degrees(longitudes) + 180.0) % 360.0 - 180.0
    columns = {}
    for index, longitude in enumerate(wrapped):
        key = round(float(longitude), 6)
        if key not in columns:
            columns[key] = index
    column_order = [columns[key] for key in sorted(columns)]
    row_order = np.argsort(90.0 - np.degrees(colatitudes))
    return values[np.ix_(row_order, column_order)], np.array(sorted(columns))


@pytest.mark.parametrize("projection, axis_name", [("mollweide", "mollweide"), ("plate_carree", "rectilinear")])
def test_matplotlib_projections_draw_the_reordered_field(projection, axis_name):
    """Matplotlib projections draw the field reordered by latitude and wrapped longitude, with a colorbar."""
    figure, axis = plot_map(LONGITUDES, COLATITUDES, FIELD, projection=projection, use_cartopy=False)
    assert axis.name == axis_name
    assert axis.figure is figure
    drawn, longitude_centers = expected_drawn_field(LONGITUDES, COLATITUDES, FIELD)
    # 0 and 2 pi are one column.
    assert drawn.shape == (7, 12)
    assert longitude_centers[0] == -180.0 and longitude_centers[-1] == 150.0
    mesh = axis.collections[0]
    np.testing.assert_allclose(np.asarray(mesh.get_array()).reshape(drawn.shape), drawn)
    assert len(figure.axes) == 2


def test_cell_edges_are_clipped_to_the_globe():
    """Cell edges are midpoints clipped to the map bounds."""
    edges = maps_module.map_cell_edges(np.array([-170.0, 0.0, 170.0]), -180.0, 180.0)
    np.testing.assert_allclose(edges, [-180.0, -85.0, 85.0, 180.0])
    np.testing.assert_allclose(maps_module.map_cell_edges(np.array([10.0]), -90.0, 90.0), [-90.0, 90.0])


# 360 repeats 0 and is dropped; on a map centered on 180 the seam sits on the map edge.
@pytest.mark.parametrize(
    "central_longitude, expected_wrapped, expected_order",
    [(0.0, [-180.0, -90.0, 0.0, 90.0], [2, 3, 0, 1]), (180.0, [0.0, 90.0, 180.0, 270.0], [0, 1, 2, 3])],
    ids=["center_0", "center_180"])
def test_longitudes_wrap_around_the_map_center(central_longitude, expected_wrapped, expected_order):
    """Longitudes wrap around the map center and a repeated 360 is dropped."""
    longitudes = np.radians([0.0, 90.0, 180.0, 270.0, 360.0])
    column_order, wrapped = maps_module.map_longitude_order(longitudes, central_longitude)
    np.testing.assert_allclose(wrapped, expected_wrapped, atol=1.0e-12)
    np.testing.assert_array_equal(column_order, expected_order)


def test_nan_cells_are_blank():
    """NaN cells are masked in the drawn mesh."""
    field = FIELD.copy()
    field[0, :] = np.nan
    _, axis = plot_map(LONGITUDES, COLATITUDES, field, use_cartopy=False)
    drawn = axis.collections[0].get_array()
    # The north-most row is drawn last.
    assert np.ma.count_masked(drawn) == 12


def test_color_scales():
    """Symmetric, explicit, and log color scales set the expected norms and masks."""
    _, axis = plot_map(LONGITUDES, COLATITUDES, FIELD, symmetric=True, use_cartopy=False)
    norm = axis.collections[0].norm
    assert -norm.vmin == norm.vmax == pytest.approx(np.max(np.abs(FIELD)))
    assert axis.collections[0].cmap.name == MAP_PLOT_STYLE["symmetric_colormap"]

    _, axis = plot_map(LONGITUDES, COLATITUDES, FIELD, value_limits=(-2.0, 3.0), use_cartopy=False)
    assert (axis.collections[0].norm.vmin, axis.collections[0].norm.vmax) == (-2.0, 3.0)

    _, axis = plot_map(LONGITUDES, COLATITUDES, FIELD, log_scale=True, colorbar=False, use_cartopy=False)
    mesh = axis.collections[0]
    assert isinstance(mesh.norm, matplotlib.colors.LogNorm)
    assert np.ma.count_masked(mesh.get_array()) == int(np.sum(expected_drawn_field(LONGITUDES, COLATITUDES, FIELD)[0]
                                                               <= 0.0))


def test_panels_from_make_map_axes():
    """plot_map draws into a given panel from make_map_axes and leaves the others empty."""
    figure, axes = make_map_axes(nrows=2, ncols=3, use_cartopy=False)
    assert axes.shape == (2, 3)
    assert all(axis.name == "mollweide" for axis in axes.flat)
    returned_figure, returned_axis = plot_map(LONGITUDES, COLATITUDES, FIELD, axis=axes[1, 2], title="Panel")
    assert returned_figure is figure and returned_axis is axes[1, 2]
    assert axes[1, 2].get_title() == "Panel"
    assert len(axes[0, 0].collections) == 0


@pytest.mark.parametrize(
    "field, kwargs",
    [
        (FIELD.T, {}),
        (FIELD, {"projection": "azimuthal"}),
        (-np.abs(FIELD), {"log_scale": True}),
        (np.full(FIELD.shape, np.nan), {}),
        (FIELD, {"symmetric": True, "log_scale": True}),
    ],
    ids=["transposed_field", "unknown_projection", "log_of_non_positive", "all_nan", "symmetric_and_log"])
def test_invalid_requests_raise(field, kwargs):
    """Invalid fields or option combinations raise ValueError."""
    with pytest.raises(ValueError):
        plot_map(LONGITUDES, COLATITUDES, field, use_cartopy=False, **kwargs)


def test_cartopy_only_requests_raise_without_cartopy():
    """Cartopy-only projections and a central longitude raise ValueError without cartopy."""
    assert MAP_PROJECTIONS["robinson"] and not MAP_PROJECTIONS["mollweide"]
    with pytest.raises(ValueError):
        make_map_axes(projection="robinson", use_cartopy=False)
    with pytest.raises(ValueError):
        make_map_axes(central_longitude=180.0, use_cartopy=False)


def test_cartopy_missing(monkeypatch):
    """Requiring missing cartopy raises ImportError; the default falls back to matplotlib."""
    monkeypatch.setattr(maps_module, "ccrs", None)
    with pytest.raises(ImportError):
        make_map_axes(use_cartopy=True)
    figure, axes = make_map_axes()
    assert axes[0, 0].name == "mollweide"


@pytest.mark.parametrize("projection", ["mollweide", "robinson", "plate_carree"])
def test_cartopy_projections(projection):
    """Cartopy projections draw every cell and render offline."""
    ccrs = pytest.importorskip("cartopy.crs")
    figure, axis = plot_map(
        LONGITUDES,
        COLATITUDES,
        FIELD,
        projection=projection,
        central_longitude=180.0,
        symmetric=True,
        use_cartopy=True)
    expected_class = {"mollweide": ccrs.Mollweide, "robinson": ccrs.Robinson, "plate_carree": ccrs.PlateCarree}
    assert isinstance(axis.projection, expected_class[projection])
    drawn, _ = expected_drawn_field(LONGITUDES, COLATITUDES, FIELD)
    mesh = axis.collections[0]
    assert np.ma.count(mesh.get_array()) == drawn.size
    # Rendering must work offline: no Natural Earth data is used.
    figure.canvas.draw()
