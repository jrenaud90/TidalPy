import matplotlib.gridspec as gspec
import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np

from TidalPy.exceptions import ParameterMissingError


def success_grid_plot(success_by_rheo: dict, x: np.ndarray, y: np.ndarray,
                      xname: str = 'unknown x', yname: str = 'unknown y',
                      xlog: bool = False, ylog: bool = False,
                      auto_show: bool = True):
    # Setup Cmap
    map_colors = ['#000000',  # Computation Fail, Code -2 --> Bad Integration
                  '#660000',  # Computation Fail, Code -1 --> Bad Year Interpolation
                  '#ffffff',  # Computation Pass, Model Fail, 0 <= Index < 1
                  '#000099',  # Computation Pass, Model Okay, 1 <= Index < 2
                  '#008000',  # Computation Pass, Model Success, Index >= 2
                  ]
    map_colors = [colors.ColorConverter.to_rgb(col) for col in map_colors]
    colors_n = len(map_colors)
    c_red = [(i / (colors_n - 1), map_colors[i][0], map_colors[i][0]) for i in range(colors_n)]
    c_green = [(i / (colors_n - 1), map_colors[i][1], map_colors[i][1]) for i in range(colors_n)]
    c_blue = [(i / (colors_n - 1), map_colors[i][2], map_colors[i][2]) for i in range(colors_n)]
    # Make the red a hard transition
    # c_red[0] = (c_red[0][0], c_red[0][1], c_red[1][1])
    # c_blue[0] = (c_blue[0][0], c_blue[0][1], c_blue[1][1])
    # c_green[0] = (c_green[0][0], c_green[0][1], c_green[1][1])
    # Make the white a hard transition
    c_red[1] = (c_red[1][0], c_red[1][1], 1.)
    c_blue[1] = (c_blue[1][0], c_blue[1][1], 1.)
    c_green[1] = (c_green[1][0], c_green[1][1], 1.)
    c_dict = {'red': tuple(c_red), 'green': tuple(c_green), 'blue': tuple(c_blue)}
    cmap = colors.LinearSegmentedColormap('twoseg_linear', c_dict)
    bounds = [-2, -1, 0, 1, 2, 4]
    norm = colors.Normalize(vmin=-2, vmax=4)

    # Setup Figure
    num_rheos = len(success_by_rheo)
    if num_rheos == 1:
        # 1 column is for the colorbar
        n_cols = 2
        width_ratios = (0.9, 0.1)
    else:
        n_cols = 3
        width_ratios = (0.45, 0.45, 0.05)
    n_rows = int(np.ceil(num_rheos / 2.))
    fig = plt.figure()
    gs = gspec.GridSpec(nrows=n_rows, ncols=n_cols, width_ratios=width_ratios, wspace=0.1, hspace=0.18)
    # Setup colorbars
    cb_axes = [fig.add_subplot(gs[i, -1]) for i in range(n_rows)]
    # Plot
    for r_i, (rheo_name, pass_index) in enumerate(success_by_rheo.items()):
        r_name = rheo_name.title()
        col_i = r_i % 2
        row_i = int(np.floor(r_i / 2.))
        ax = fig.add_subplot(gs[row_i, col_i])

        if ylog and not xlog:
            dx = (x[1] - x[0]) / 2.
            dy = np.sqrt(y[1] / y[0])
            ax.set_yscale('log')
            extent_x = np.linspace(x[0] - dx, x[-1] + dx, len(x) + 1)
            extent_y = np.logspace(np.log10(y[0] / dy), np.log10(y[-1] * dy), len(y) + 1)
        elif not ylog and xlog:
            dx = np.sqrt(x[1] / x[0])
            dy = (y[1] - y[0]) / 2.
            ax.set_xscale('log')
            extent_y = np.linspace(y[0] - dy, y[-1] + dy, len(y) + 1)
            extent_x = np.logspace(np.log10(x[0] / dx), np.log10(x[-1] * dx), len(x) + 1)
        elif ylog and xlog:
            dx = np.sqrt(x[1] / x[0])
            dy = np.sqrt(y[1] / y[0])
            ax.set_xscale('log')
            ax.set_yscale('log')
            extent_x = np.logspace(np.log10(x[0] / dx), np.log10(x[-1] * dx), len(x) + 1)
            extent_y = np.logspace(np.log10(y[0] / dy), np.log10(y[-1] * dy), len(y) + 1)
        else:
            dx = (x[1] - x[0]) / 2.
            dy = (y[1] - y[0]) / 2.
            extent_x = np.linspace(x[0] - dx, x[-1] + dx, len(x) + 1)
            extent_y = np.linspace(y[0] - dy, y[-1] + dy, len(y) + 1)

        # This extention of pass fail (and the X, Y len(X or Y) + 1) correct the issue of pcolormesh() not plotting the last column and row of C.
        pass_index = np.concatenate((pass_index, np.zeros((pass_index.shape[1], 1))), axis=1)
        pass_index = np.concatenate((pass_index, np.zeros((1, pass_index.shape[1]))))

        cb = ax.pcolormesh(extent_x, extent_y, pass_index, cmap=cmap, norm=norm)
        ax.set_title(r_name)
        if col_i == 0:
            ax.set_ylabel(yname)
        else:
            ax.set_yticklabels([])
        if row_i == n_rows - 1:
            ax.set_xlabel(xname)
        else:
            ax.set_xticklabels([])
        ax.set_xlim((extent_x[0], extent_x[-1]))
        ax.set_ylim((extent_y[0], extent_y[-1]))

        for _x in x[1:]:
            if not xlog:
                ax.axvline(x=_x - dx, c='k', ls='-', alpha=0.7, linewidth=0.5)
            else:
                ax.axvline(x=_x / dx, c='k', ls='-', alpha=0.7, linewidth=0.5)

        for _y in y[1:]:
            if not ylog:
                ax.axhline(y=_y - dy, c='k', ls='-', alpha=0.7, linewidth=0.5)
            else:
                ax.axhline(y=_y / dy, c='k', ls='-', alpha=0.7, linewidth=0.5)

    for cb_ax in cb_axes:
        cb_obj = plt.colorbar(cb, cax=cb_ax)
        cb_obj.set_label('Pass Index')

    if auto_show:
        plt.show()

    return fig