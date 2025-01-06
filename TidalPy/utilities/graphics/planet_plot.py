import numpy as np
from matplotlib import gridspec as gspec, pyplot as plt

from TidalPy.exceptions import ParameterMissingError

SCATTER_SIZE = 35
SCATTER_SHAPE_1 = "."
SCATTER_SHAPE_2 = "x"
LS_1 = '-'
LS_2 = ':'
SCALE = 4

GRAVITY_COLOR     = 'g'
DENSITY_COLOR     = 'k'
PRESSURE_COLOR    = 'b'
TEMPERATURE_COLOR = 'orange'
SHEAR_COLOR       = 'm'
BULK_COLOR        = 'r'

FONTSIZE_1 = 12
FONTSIZE_2 = 14

def planet_plot(
    radii: np.ndarray,
    gravity_array: np.ndarray,
    pressure_array: np.ndarray,
    density_array: np.ndarray,
    temperature_array: np.ndarray = None,
    shear_modulus_array: np.ndarray = None,
    bulk_modulus_array: np.ndarray = None,
    planet_radius: float = None,
    bulk_density: float = None,
    planet_name: str = None,
    use_scatter: bool = False,
    depth_plot: bool = False,
    auto_show: bool = False,
    annotate: bool = True
    ):
    """ Plots the depth plot of a planet in 3 or 4 panels (temperature is optional)

    Parameters
    ----------
    radii : np.ndarray
        Planet's radius discretized into multiple chunks in a 1-D array [m]
    gravity_array : np.ndarray
        Planet's gravity as a function of depth or radius
    pressure_array : np.ndarray
        Planet's pressure as a function of depth or radius
    density_array : np.ndarray
        Planet's density as a function of depth or radius
    temperature_array : np.ndarray
        (optional) Planet's internal temperature as a function of depth or radius
    shear_modulus_array : np.ndarray
        (optional) Shear modulus as a function of depth or radius
        Can be complex-valued.
    bulk_modulus_array : np.ndarray
        (optional) Shear modulus as a function of depth or radius
        Can be complex-valued.
    planet_radius : float
        (optional) Planet's outer radius - used for depth plot
    bulk_density : float
        (optional) Planet's bulk density can be shown as an annotation
    planet_name : str
        Planet's name
    use_scatter : bool
        (default = False) If True, then scatter rather than line plot will be shown.
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
        y = planet_radius - radii
    else:
        y = radii
    shape = y.shape

    assert gravity_array.shape == shape
    assert pressure_array.shape == shape
    assert density_array.shape == shape
    use_temperature = False
    if temperature_array is not None:
        assert temperature_array.shape == shape
        use_temperature = True
    use_modulus = False
    if shear_modulus_array is not None or bulk_modulus_array is not None:
        use_modulus = True
    use_complex_modulus = False
    two_modulus = False
    if shear_modulus_array is not None:
        if shear_modulus_array.dtype == np.complex128:
            use_complex_modulus = True
    if bulk_modulus_array is not None:
        if bulk_modulus_array.dtype == np.complex128:
            use_complex_modulus = True
    if shear_modulus_array is not None and bulk_modulus_array is not None:
        two_modulus = True

    # Perform any unit conversions
    pressure_array = pressure_array / 1.e9   # Use GPa
    density_array = density_array / 1000.    # Use kg m-3
    y /= 1000.                            # Use km
    if bulk_density is not None:
        bulk_density = bulk_density / 1000. # Use kg m-3
    if shear_modulus_array is not None:
        shear_modulus_array = shear_modulus_array / 1.0e9
    if bulk_modulus_array is not None:
        bulk_modulus_array = bulk_modulus_array / 1.0e9

    n = 2
    if use_modulus:
        n = 3

    # Construct plot grid.
    fig, axes = plt.subplots(figsize=(n * SCALE, 1 * SCALE), nrows=1, ncols=n)
    fig.subplots_adjust(wspace=0.2, hspace=0.1)

    # Subplot 1
    ax_gravity = axes[0]
    ax_density = ax_gravity.twiny()
    if use_scatter:
        ax_gravity.scatter(gravity_array, y, s=SCATTER_SIZE, c=GRAVITY_COLOR, marker=SCATTER_SHAPE_1)
        ax_density.scatter(density_array, y, s=SCATTER_SIZE, c=DENSITY_COLOR, marker=SCATTER_SHAPE_1)
    else:
        ax_gravity.plot(gravity_array, y, ls=LS_1, c=GRAVITY_COLOR)
        ax_density.plot(density_array, y, ls=LS_1, c=DENSITY_COLOR)

    # Subplot 2
    ax_pressure = axes[1]
    ax_temperature = None
    if use_temperature:
        ax_temperature = ax_pressure.twiny()
    if use_scatter:
        ax_pressure.scatter(pressure_array, y, s=SCATTER_SIZE, c=PRESSURE_COLOR, marker=SCATTER_SHAPE_1)
        if ax_temperature is not None:
            ax_temperature.scatter(temperature_array, y, s=SCATTER_SIZE, c=TEMPERATURE_COLOR, marker=SCATTER_SHAPE_1)
    else:
        ax_pressure.plot(pressure_array, y, ls=LS_1, c=PRESSURE_COLOR)
        if ax_temperature is not None:
            ax_temperature.plot(temperature_array, y, ls=LS_1, c=TEMPERATURE_COLOR)
    
    # Subplot 3
    ax_modulus = None
    ax_imag_modulus = None
    if use_modulus:
        ax_modulus = axes[2]
        if use_complex_modulus:
            ax_imag_modulus = ax_modulus.twiny()
    if use_scatter:
        if shear_modulus_array is not None:
            if shear_modulus_array.dtype == np.complex128:
                ax_modulus.scatter(shear_modulus_array.real, y, s=SCATTER_SIZE, c=SHEAR_COLOR, marker=SCATTER_SHAPE_1, label="Shear")
                ax_imag_modulus.scatter(shear_modulus_array.imag, y, s=SCATTER_SIZE, c=SHEAR_COLOR, marker=SCATTER_SHAPE_2, label="Shear")
            else:
                ax_modulus.scatter(shear_modulus_array, y, s=SCATTER_SIZE, c=SHEAR_COLOR, marker=SCATTER_SHAPE_1, label="Shear")
        if bulk_modulus_array is not None:
            if bulk_modulus_array.dtype == np.complex128:
                ax_modulus.scatter(bulk_modulus_array.real, y, s=SCATTER_SIZE, c=BULK_COLOR, marker=SCATTER_SHAPE_1, label="Bulk")
                ax_imag_modulus.scatter(bulk_modulus_array.imag, y, s=SCATTER_SIZE, c=BULK_COLOR, marker=SCATTER_SHAPE_2, label="Bulk")
            else:
                ax_modulus.scatter(bulk_modulus_array, y, s=SCATTER_SIZE, c=BULK_COLOR, marker=SCATTER_SHAPE_1, label="Bulk")
    else:
        if shear_modulus_array is not None:
            if shear_modulus_array.dtype == np.complex128:
                ax_modulus.plot(shear_modulus_array.real, y, ls=LS_1, c=SHEAR_COLOR, label="Shear")
                ax_imag_modulus.plot(shear_modulus_array.imag, y, ls=LS_2, c=SHEAR_COLOR, label="Shear")
            else:
                ax_modulus.plot(shear_modulus_array, y, ls=LS_1, c=SHEAR_COLOR, label="Shear")
        if bulk_modulus_array is not None:
            if bulk_modulus_array.dtype == np.complex128:
                ax_modulus.plot(bulk_modulus_array.real, y, ls=LS_1, c=BULK_COLOR, label="Bulk")
                ax_imag_modulus.plot(bulk_modulus_array.imag, y, ls=LS_2, c=BULK_COLOR, label="Bulk")
            else:
                ax_modulus.plot(bulk_modulus_array, y, ls=LS_1, c=BULK_COLOR, label="Bulk")

    # Setup labels
    if depth_plot:
        ax_gravity.set_ylabel('Depth [km]', fontsize=FONTSIZE_1)
    else:
        ax_gravity.set_ylabel('Radius [km]', fontsize=FONTSIZE_1)
    
    ax_gravity.set_xlabel('Gravity [m s$^{-2}$]', color=GRAVITY_COLOR, fontsize=FONTSIZE_1)
    ax_density.set_xlabel('Density [kg m$^{-3}$]', color=DENSITY_COLOR, fontsize=FONTSIZE_1)
    ax_pressure.set_xlabel('Pressure [GPa]', color=PRESSURE_COLOR, fontsize=FONTSIZE_1)
    if ax_temperature is not None:
        ax_temperature.set_xlabel('Temperature [K]', color=TEMPERATURE_COLOR, fontsize=FONTSIZE_1)
    if ax_modulus is not None:
        if use_complex_modulus:
            if use_scatter:
                ax_modulus.set_xlabel('Re[Modulus] (Dots) [GPa]', color='k', fontsize=FONTSIZE_1)
            else:
                ax_modulus.set_xlabel('Re[Modulus] (Solid) [GPa]', color='k', fontsize=FONTSIZE_1)
        else:
            ax_modulus.set_xlabel('Modulus [GPa]', color='k', fontsize=FONTSIZE_1)
    if ax_imag_modulus is not None:
        if use_scatter:
            ax_imag_modulus.set_xlabel('Im[Modulus] (Xs) [GPa]', color='k', fontsize=FONTSIZE_1)
        else:
            ax_imag_modulus.set_xlabel('Im[Modulus] (Dotted) [GPa]', color='k', fontsize=FONTSIZE_1)
        ax_imag_modulus.set_xscale('log')

    # Add annotations
    if annotate:
        if depth_plot:
            grav_pos = (0.05, 0.9)
        else:
            grav_pos = (0.25, 0.05)
        ax_gravity.text(
            *grav_pos,
            '$g_{s}$' + f' = {gravity_array[-1]:0.2f}',
            horizontalalignment='left', verticalalignment='center', transform=ax_gravity.transAxes
        )

        if depth_plot:
            press_pos = (0.05, 0.15)
        else:
            press_pos = (0.05, 0.05)
        ax_pressure.text(
            *press_pos,
            '$P_{0}$' + f' = {pressure_array[0]:0.2f}',
            horizontalalignment='left', verticalalignment='center', transform=ax_pressure.transAxes
        )

        if bulk_density is not None:
            if depth_plot:
                density_pos = (0.05, 0.15)
            else:
                density_pos = (0.25, 0.90)
            ax_density.text(
                *density_pos,
                '$\\bar{\\rho}$' + f' = {bulk_density:0.2f}',
                horizontalalignment='left', verticalalignment='center', transform=ax_density.transAxes
            )
    
    # Setup Legend
    if two_modulus:
        ax_modulus.legend(loc='best')
    
    # Remove uneeded tick labels
    ax_pressure.set_yticklabels([])
    if ax_temperature is not None:
        ax_temperature.set_yticklabels([])
    if ax_modulus is not None:
        ax_modulus.set_yticklabels([])
    if ax_imag_modulus is not None:
        ax_imag_modulus.set_yticklabels([])

    # Other parameters
    if planet_name is not None:
        fig.suptitle(planet_name.title(), y=1.1, fontsize=FONTSIZE_2)

    if auto_show:
        plt.show()

    return fig, axes
