"""Tests for the surface map helpers in TidalPy.Utilities_x.graphics_x (`plot_map`, `make_map_axes`).

Figures are drawn on the non-interactive Agg backend and closed after each test. The matplotlib projections are
tested everywhere; the cartopy projections only where cartopy is installed.
"""
import matplotlib

matplotlib.use("Agg")

import numpy as np
import pytest
from matplotlib import pyplot as plt

from TidalPy.Utilities_x.graphics_x import MAP_PLOT_STYLE, MAP_PROJECTIONS, make_map_axes, plot_map
from TidalPy.Utilities_x.graphics_x import maps as maps_module

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


# =====================================================================================================================
# Matplotlib projections
# =====================================================================================================================
@pytest.mark.parametrize("projection, axis_name", [("mollweide", "mollweide"), ("plate_carree", "rectilinear")])
def test_matplotlib_projections_draw_the_reordered_field(projection, axis_name):
    figure, axis = plot_map(LONGITUDES, COLATITUDES, FIELD, projection=projection, use_cartopy=False)
    assert axis.name == axis_name
    assert axis.figure is figure
    drawn, longitude_centers = expected_drawn_field(LONGITUDES, COLATITUDES, FIELD)
    assert drawn.shape == (7, 12)                       # 0 and 2 pi are one column
    assert longitude_centers[0] == -180.0 and longitude_centers[-1] == 150.0
    mesh = axis.collections[0]
    np.testing.assert_allclose(np.asarray(mesh.get_array()).reshape(drawn.shape), drawn)
    assert len(figure.axes) == 2                        # The panel and its colorbar


def test_cell_edges_are_clipped_to_the_globe():
    edges = maps_module.map_cell_edges(np.array([-170.0, 0.0, 170.0]), -180.0, 180.0)
    np.testing.assert_allclose(edges, [-180.0, -85.0, 85.0, 180.0])
    np.testing.assert_allclose(maps_module.map_cell_edges(np.array([10.0]), -90.0, 90.0), [-90.0, 90.0])


def test_longitudes_wrap_around_the_map_center():
    longitudes = np.radians([0.0, 90.0, 180.0, 270.0, 360.0])
    column_order, wrapped = maps_module.map_longitude_order(longitudes, 0.0)
    np.testing.assert_allclose(wrapped, [-180.0, -90.0, 0.0, 90.0], atol=1.0e-12)
    np.testing.assert_array_equal(column_order, [2, 3, 0, 1])       # 360 repeats 0 and is dropped
    column_order, wrapped = maps_module.map_longitude_order(longitudes, 180.0)
    np.testing.assert_allclose(wrapped, [0.0, 90.0, 180.0, 270.0], atol=1.0e-12)
    np.testing.assert_array_equal(column_order, [0, 1, 2, 3])       # The seam sits on the edge of a map centered on 180


def test_nan_cells_are_blank():
    field = FIELD.copy()
    field[0, :] = np.nan
    _, axis = plot_map(LONGITUDES, COLATITUDES, field, use_cartopy=False)
    drawn = axis.collections[0].get_array()
    assert np.ma.count_masked(drawn) == 12              # The north-most row, drawn last, is fully masked


def test_color_scales():
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
    figure, axes = make_map_axes(nrows=2, ncols=3, use_cartopy=False)
    assert axes.shape == (2, 3)
    assert all(axis.name == "mollweide" for axis in axes.flat)
    returned_figure, returned_axis = plot_map(LONGITUDES, COLATITUDES, FIELD, axis=axes[1, 2], title="Panel")
    assert returned_figure is figure and returned_axis is axes[1, 2]
    assert axes[1, 2].get_title() == "Panel"
    assert len(axes[0, 0].collections) == 0


# =====================================================================================================================
# Validation and backend selection
# =====================================================================================================================
def test_invalid_requests_raise():
    with pytest.raises(ValueError):
        plot_map(LONGITUDES, COLATITUDES, FIELD.T, use_cartopy=False)
    with pytest.raises(ValueError):
        plot_map(LONGITUDES, COLATITUDES, FIELD, projection="azimuthal", use_cartopy=False)
    with pytest.raises(ValueError):
        plot_map(LONGITUDES, COLATITUDES, -np.abs(FIELD), log_scale=True, use_cartopy=False)
    with pytest.raises(ValueError):
        plot_map(LONGITUDES, COLATITUDES, np.full(FIELD.shape, np.nan), use_cartopy=False)
    with pytest.raises(ValueError):
        plot_map(LONGITUDES, COLATITUDES, FIELD, symmetric=True, log_scale=True, use_cartopy=False)


def test_cartopy_only_requests_raise_without_cartopy():
    assert MAP_PROJECTIONS["robinson"] and not MAP_PROJECTIONS["mollweide"]
    with pytest.raises(ValueError):
        make_map_axes(projection="robinson", use_cartopy=False)
    with pytest.raises(ValueError):
        make_map_axes(central_longitude=180.0, use_cartopy=False)


def test_cartopy_missing(monkeypatch):
    monkeypatch.setattr(maps_module, "ccrs", None)
    with pytest.raises(ImportError):
        make_map_axes(use_cartopy=True)
    figure, axes = make_map_axes()                      # The default falls back to matplotlib
    assert axes[0, 0].name == "mollweide"


# =====================================================================================================================
# Cartopy projections
# =====================================================================================================================
@pytest.mark.parametrize("projection", ["mollweide", "robinson", "plate_carree"])
def test_cartopy_projections(projection):
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
    figure.canvas.draw()                                # Rendering must work offline: no Natural Earth data is used
