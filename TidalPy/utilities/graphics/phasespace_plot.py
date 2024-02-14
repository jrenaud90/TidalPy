import numpy as np
from matplotlib import colors, gridspec as gspec, pyplot as plt

from TidalPy.utilities.numpy_helper.array_other import normalize_dict


def success_grid_plot(
    success_by_rheo: dict, x: np.ndarray, y: np.ndarray,
    xname: str = 'unknown x', yname: str = 'unknown y',
    xlog: bool = False, ylog: bool = False,
    print_value: bool = False,
    auto_show: bool = True
    ):
    cb_N = 7
    # Normalize input data
    success_by_rheo_norm = normalize_dict(success_by_rheo, pass_negatives=True, new_max=4.0, new_min=0.0)
    new_success_by_rheo = dict()
    for rheo, data in success_by_rheo_norm.items():
        # Make the fail codes fall inbetween the ranges to ensure the get the correct color below
        new_data = data.copy()
        new_data[new_data == -2] = -1.5
        new_data[new_data == -1] = -0.5
        new_success_by_rheo[rheo] = new_data

    # Setup Cmap
    map_colors = ['#000000',  # Computation Fail, Code -2 --> Bad Integration
                  '#660000',  # Computation Fail, Code -1 --> Bad Year Interpolation
                  '#ffffff',  # Computation Pass, Model Fail, 0 <= Index < 1
                  '#000099',  # Computation Pass, Model Okay, 1 <= Index < 2
                  '#008000',  # Computation Pass, Model Success, Index >= 2
                  ]
    map_colors = [colors.ColorConverter.to_rgb(col) for col in map_colors]
    c_dict = {
        'red'  : ((0, map_colors[0][0], map_colors[0][0]),  # Start of -2 Code
                  (1 / (cb_N - 1), map_colors[0][0], map_colors[0][0]),  # Hard Stop of -2 Code
                  (1 / (cb_N - 1), map_colors[1][0], map_colors[1][0]),  # Hard Start of -1 Code
                  (1.99 / (cb_N - 1), map_colors[1][0], map_colors[1][0]),  # Hard Stop of -1 Code
                  (1.99 / (cb_N - 1), map_colors[3][0], map_colors[2][0]),
                  # These make no sense to me, but trial and error found them...
                  (6 / (cb_N - 1), map_colors[4][0], map_colors[3][0])
                  ),
        'green': ((0, map_colors[0][1], map_colors[0][1]),  # Start of -2 Code
                  (1 / (cb_N - 1), map_colors[0][1], map_colors[0][1]),  # Hard Stop of -2 Code
                  (1 / (cb_N - 1), map_colors[1][1], map_colors[1][1]),  # Hard Start of -1 Code
                  (1.99 / (cb_N - 1), map_colors[1][1], map_colors[1][1]),  # Hard Stop of -1 Code
                  (1.99 / (cb_N - 1), map_colors[3][1], map_colors[2][1]),
                  (6 / (cb_N - 1), map_colors[4][1], map_colors[3][1])
                  ),
        'blue' : ((0, map_colors[0][2], map_colors[0][2]),  # Start of -2 Code
                  (1 / (cb_N - 1), map_colors[0][2], map_colors[0][2]),  # Hard Stop of -2 Code
                  (1 / (cb_N - 1), map_colors[1][2], map_colors[1][2]),  # Hard Start of -1 Code
                  (1.99 / (cb_N - 1), map_colors[1][2], map_colors[1][2]),  # Hard Stop of -1 Code
                  (1.99 / (cb_N - 1), map_colors[3][2], map_colors[2][2]),
                  (6 / (cb_N - 1), map_colors[4][2], map_colors[3][2])
                  ),
        }

    cmap = colors.LinearSegmentedColormap('twoseg_linear', c_dict)
    norm = colors.Normalize(vmin=-2, vmax=4)

    # Setup Figure
    num_rheos = len(new_success_by_rheo)
    if num_rheos == 1:
        # 1 column is for the colorbar
        n_cols = 2
        width_ratios = (0.9, 0.1)
    else:
        n_cols = 3
        width_ratios = (0.45, 0.45, 0.05)
    n_rows = int(np.ceil(num_rheos / 2.))
    fig = plt.figure()
    gs = gspec.GridSpec(nrows=n_rows, ncols=n_cols, width_ratios=width_ratios, wspace=0.1, hspace=0.22)
    # Setup colorbars
    cb_axes = [fig.add_subplot(gs[i, -1]) for i in range(n_rows)]
    # Plot
    for r_i, (rheo_name, pass_index) in enumerate(new_success_by_rheo.items()):
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

        # TODO: I don't think these concats are actually needed
        # This extension of pass_index (and the X, Y len(X or Y) + 1) corrects the issue of pcolormesh() not plotting the last column and row of C.
        # pass_index = np.concatenate((pass_index, np.zeros((pass_index.shape[1], 1))), axis=1)
        # pass_index = np.concatenate((pass_index, np.zeros((1, pass_index.shape[1]))))

        # The edgecolor and linewidth control the gridlines
        cb = ax.pcolormesh(
            extent_x, extent_y, np.transpose(pass_index), cmap=cmap, norm=norm, edgecolor='k',
            linewidths=.5
            )
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

        # Manually set x and y ticks to match the inputs since each box corresponds to a specific case
        ax.xaxis.set_ticks(x)
        ax.yaxis.set_ticks(y)

        if print_value:
            # Print the success value on the middle of the box.
            x_offset = 0.1 * dx
            y_offset = 0.2 * dy
            for i in range(len(extent_x) - 1):
                for j in range(len(extent_y) - 1):
                    x_, y_ = extent_x[i], extent_y[j]
                    val = success_by_rheo[rheo_name][i, j]
                    corrected_vale = new_success_by_rheo[rheo_name][i, j]
                    text = f'{val:0.03f}, {corrected_vale:0.03f}'
                    if val < 0:
                        col = 'white'
                    else:
                        col = 'black'

                    ax.text(x_ + x_offset, y_ + y_offset, text, size=8, alpha=0.7, color=col)

    for cb_ax in cb_axes:
        cb_obj = plt.colorbar(cb, cax=cb_ax, ticks=[-1.75, -0.75, 0, 4])
        cb_obj.set_label('Pass Index')
        cb_ax.set_yticklabels(['E1', 'E2', '$0$', '$1$'])

    if auto_show:
        plt.show()

    return fig
