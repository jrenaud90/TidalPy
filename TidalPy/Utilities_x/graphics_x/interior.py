"""Plot a planet's interior profiles as found by the equation-of-state solver.

`plot_interior` draws gravity and density (shared panel, two x-axes), pressure with optional temperature,
and, when given, the shear and bulk moduli (real parts, plus imaginary parts on a twin axis for complex
moduli) against radius or depth.
"""

from __future__ import annotations

from typing import Dict, Optional, Tuple

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.figure import Figure

# =====================================================================================================================
# Style
# =====================================================================================================================

# Colors, line styles, markers, and sizes used by `plot_interior`; edit in place to restyle every plot.
INTERIOR_PLOT_STYLE: Dict[str, object] = {
    "gravity_color": "g",
    "density_color": "k",
    "pressure_color": "b",
    "temperature_color": "orange",
    "shear_color": "m",
    "bulk_color": "r",
    "line_style": "-",
    "imaginary_line_style": ":",
    "marker": ".",
    "imaginary_marker": "x",
    "marker_size": 35,
    "panel_size_inches": 4.0,
    "label_fontsize": 12,
    "title_fontsize": 14,
}


# =====================================================================================================================
# Plotting
# =====================================================================================================================

def _draw(axis, values, vertical, color, use_scatter: bool, imaginary: bool = False, label: Optional[str] = None):
    """Draw one profile as a line or scatter series."""
    style = INTERIOR_PLOT_STYLE
    if use_scatter:
        axis.scatter(values, vertical, s=style["marker_size"], c=color,
                     marker=style["imaginary_marker"] if imaginary else style["marker"], label=label)
    else:
        axis.plot(values, vertical, ls=style["imaginary_line_style"] if imaginary else style["line_style"], c=color,
                  label=label)


def plot_interior(
        radius: np.ndarray,
        gravity: np.ndarray,
        pressure: np.ndarray,
        density: np.ndarray,
        temperature: Optional[np.ndarray] = None,
        shear_modulus: Optional[np.ndarray] = None,
        bulk_modulus: Optional[np.ndarray] = None,
        planet_radius: Optional[float] = None,
        bulk_density: Optional[float] = None,
        planet_name: Optional[str] = None,
        depth_plot: bool = False,
        use_scatter: bool = False,
        annotate: bool = True,
        show_plot: bool = False,
        ) -> Tuple[Figure, np.ndarray]:
    """Plot a planet's interior profiles in two or three panels.

    Parameters
    ----------
    radius : array
        Radius of each point [m].
    gravity : array
        Gravitational acceleration at each radius [m s-2].
    pressure : array
        Pressure at each radius [Pa]; plotted in GPa.
    density : array
        Density at each radius [kg m-3].
    temperature : array, optional
        Temperature at each radius [K]; drawn on a twin axis of the pressure panel.
    shear_modulus, bulk_modulus : array, optional
        Shear and bulk moduli at each radius [Pa], real or complex; plotted in GPa in a third panel. The
        imaginary parts of complex moduli go on a twin axis (log scale when they are all positive).
    planet_radius : float, optional
        Planet radius [m]; required when `depth_plot` is True.
    bulk_density : float, optional
        Planet bulk density [kg m-3], shown as an annotation.
    planet_name : str, optional
        Figure title.
    depth_plot : bool, default False
        Plot against depth (``planet_radius - radius``) instead of radius.
    use_scatter : bool, default False
        Draw points instead of lines.
    annotate : bool, default True
        Annotate the surface gravity, central pressure, and bulk density.
    show_plot : bool, default False
        Call ``matplotlib.pyplot.show()`` before returning.

    Returns
    -------
    figure : matplotlib.figure.Figure
    axes : numpy.ndarray of matplotlib.axes.Axes
        The primary panels: gravity/density, pressure(/temperature), and moduli when given.

    Assumptions
    -----------
    - Arrays are ordered from the planet center outward, so the last gravity value is the surface gravity
      and the first pressure value is the central pressure.
    """
    style = INTERIOR_PLOT_STYLE
    radius = np.asarray(radius, dtype=np.float64).ravel()
    if depth_plot and planet_radius is None:
        raise ValueError("`planet_radius` is required for a depth plot.")
    profiles = {"gravity": gravity, "pressure": pressure, "density": density}
    if temperature is not None:
        profiles["temperature"] = temperature
    if shear_modulus is not None:
        profiles["shear_modulus"] = shear_modulus
    if bulk_modulus is not None:
        profiles["bulk_modulus"] = bulk_modulus
    arrays = {}
    for name, values in profiles.items():
        values = np.asarray(values).ravel()
        if values.size != radius.size:
            raise ValueError(f"`{name}` has {values.size} points but `radius` has {radius.size}.")
        arrays[name] = values

    vertical = (planet_radius - radius) if depth_plot else radius
    vertical_km = vertical / 1000.0
    pressure_gpa = arrays["pressure"] / 1.0e9
    use_moduli = shear_modulus is not None or bulk_modulus is not None
    moduli = {name: arrays[name] / 1.0e9 for name in ("shear_modulus", "bulk_modulus") if name in arrays}
    use_complex = any(np.iscomplexobj(values) for values in moduli.values())

    panel_count = 3 if use_moduli else 2
    size = float(style["panel_size_inches"])
    figure, axes = plt.subplots(nrows=1, ncols=panel_count, figsize=(panel_count * size, size))
    figure.subplots_adjust(wspace=0.2, hspace=0.1)

    # Panel 1: gravity with density on a twin axis.
    ax_gravity = axes[0]
    ax_density = ax_gravity.twiny()
    _draw(ax_gravity, arrays["gravity"], vertical_km, style["gravity_color"], use_scatter)
    _draw(ax_density, arrays["density"], vertical_km, style["density_color"], use_scatter)
    ax_gravity.set_ylabel("Depth [km]" if depth_plot else "Radius [km]", fontsize=style["label_fontsize"])
    ax_gravity.set_xlabel("Gravity [m s$^{-2}$]", color=style["gravity_color"], fontsize=style["label_fontsize"])
    ax_density.set_xlabel("Density [kg m$^{-3}$]", color=style["density_color"], fontsize=style["label_fontsize"])

    # Panel 2: pressure with optional temperature.
    ax_pressure = axes[1]
    _draw(ax_pressure, pressure_gpa, vertical_km, style["pressure_color"], use_scatter)
    ax_pressure.set_xlabel("Pressure [GPa]", color=style["pressure_color"], fontsize=style["label_fontsize"])
    ax_pressure.yaxis.set_ticklabels([])
    if "temperature" in arrays:
        ax_temperature = ax_pressure.twiny()
        _draw(ax_temperature, arrays["temperature"], vertical_km, style["temperature_color"], use_scatter)
        ax_temperature.set_xlabel("Temperature [K]", color=style["temperature_color"],
                                  fontsize=style["label_fontsize"])

    # Panel 3: moduli (real parts, imaginary parts on a twin axis).
    if use_moduli:
        ax_modulus = axes[2]
        ax_imaginary = ax_modulus.twiny() if use_complex else None
        colors = {"shear_modulus": style["shear_color"], "bulk_modulus": style["bulk_color"]}
        names = {"shear_modulus": "Shear", "bulk_modulus": "Bulk"}
        imaginary_values = []
        for name, values in moduli.items():
            _draw(ax_modulus, np.real(values), vertical_km, colors[name], use_scatter, label=names[name])
            if ax_imaginary is not None and np.iscomplexobj(values):
                _draw(ax_imaginary, np.imag(values), vertical_km, colors[name], use_scatter, imaginary=True)
                imaginary_values.append(np.imag(values))
        marker_note = "points" if use_scatter else "solid"
        ax_modulus.set_xlabel(f"Re[Modulus] ({marker_note}) [GPa]" if use_complex else "Modulus [GPa]",
                              fontsize=style["label_fontsize"])
        ax_modulus.yaxis.set_ticklabels([])
        if len(moduli) > 1:
            ax_modulus.legend(loc="best")
        if ax_imaginary is not None:
            imaginary_note = "crosses" if use_scatter else "dotted"
            ax_imaginary.set_xlabel(f"Im[Modulus] ({imaginary_note}) [GPa]", fontsize=style["label_fontsize"])
            stacked = np.concatenate(imaginary_values)
            if stacked.size > 0 and np.all(stacked > 0.0):
                ax_imaginary.set_xscale("log")

    # Annotations.
    if annotate:
        ax_gravity.text(*((0.05, 0.90) if depth_plot else (0.25, 0.05)),
                        "$g_{s}$" + f" = {arrays['gravity'][-1]:0.2f} m s$^{{-2}}$",
                        horizontalalignment="left", verticalalignment="center", transform=ax_gravity.transAxes)
        ax_pressure.text(*((0.05, 0.15) if depth_plot else (0.05, 0.05)),
                         "$P_{0}$" + f" = {pressure_gpa[0]:0.2f} GPa",
                         horizontalalignment="left", verticalalignment="center", transform=ax_pressure.transAxes)
        if bulk_density is not None:
            ax_density.text(*((0.05, 0.15) if depth_plot else (0.25, 0.90)),
                            "$\\bar{\\rho}$" + f" = {bulk_density:0.1f} kg m$^{{-3}}$",
                            horizontalalignment="left", verticalalignment="center", transform=ax_density.transAxes)

    if planet_name is not None:
        figure.suptitle(planet_name.title(), y=1.1, fontsize=style["title_fontsize"])
    if show_plot:
        plt.show()
    return figure, axes
