import importlib.util
from typing import Dict, List, Tuple, Union

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from TidalPy.paths import unique_path
from TidalPy.utilities.classes.base import TidalPyClass

# Cartopy is usually not installed during testing which is fine. But this module will throw an error when tests are run
#    and Cartopy is not installed. So check if it is installed before using projections.
CARTOPY_INSTALLED = False
KNOWN_PROJECTIONS = None
if importlib.util.find_spec('cartopy') is not None:
    spec = importlib.util.find_spec('cartopy')
    from .global_map import KNOWN_PROJECTIONS
    CARTOPY_INSTALLED = True

AxisReferencedInputStr = Union[str,
                               List[str],
                               Tuple[str, ...],
                               Dict[Tuple[int, int], str]]
AxisReferencedInputBool = Union[bool,
                                List[bool],
                                Tuple[bool, ...],
                                Dict[Tuple[int, int], bool]]
AxisReferencedInputLimit = Union[Tuple[float, float],
                                 List[Tuple[float, float]],
                                 Tuple[Tuple[float, float], ...],
                                 Dict[Tuple[int, int], Tuple[float, float]]]
AxisReferencedInput = Union[AxisReferencedInputStr, AxisReferencedInputBool, AxisReferencedInputLimit]
AxisReferencedOutput = Union[Dict[Tuple[int, int], str],
                             Dict[Tuple[int, int], bool],
                             Dict[Tuple[int, int], Tuple[float, float]]]


class GridPlot(TidalPyClass):

    def __init__(
        self, nrows: int, ncols: int,
        make_colorbars: Union[bool, List[bool]] = False, colorbar_ratio: float = 0.1,
        projections: Union[str, List[str]] = None,
        figure_scale: float = 1., aspect_ratio: float = 2.,
        hspace: float = 0.25, wspace: float = 0.3, use_tight_layout: bool = True,
        suptitle: str = None, titles: AxisReferencedInputStr = None,
        xlabels: AxisReferencedInputStr = None, xscales: AxisReferencedInputStr = None,
        ylabels: AxisReferencedInputStr = None, yscales: AxisReferencedInputStr = None,
        label_kwargs: dict = None, dark_mode: bool = False
        ):
        """ Helper function to quickly make a figure of plots in a grid.

        Built off of matplotlib's gridspec

        Parameters
        ----------
        nrows : int
            Number of rows in the grid
        ncols : int
            Number of columns in the grid
        make_colorbars : bool = False
            If `True`, then an additional set of colorbar columns will be added to the plots
        colorbar_ratio : float = 0.1
            Relative size of the colorbar axes compared to the figure size
        projections : Union[str, List[str]] = None
            List of optional projections applied to each subplot
        figure_scale : float = 1.
            Size of the figure
        aspect_ratio : float = 2.
            Aspect ratio of the figure
        hspace : float = 0.25
            Vertical space between subplots
        wspace : float = 0.3
            Horizontal space between subplots
        use_tight_layout : bool = True
            If `True`, then matplotlib's tight_layout feature will be used
        suptitle : str = None
            Optional super title for the figure
        titles : AxisReferencedInputStr = None
            Optional list of titles for each subplot
        xlabels : AxisReferencedInputStr = None
            Optional list of x labels for each subplot
        xscales : AxisReferencedInputStr = None
            Optional list of x scales for each subplot
        ylabels : AxisReferencedInputStr = None
            Optional list of y labels for each subplot
        yscales : AxisReferencedInputStr = None
            Optional list of y scales for each subplot
        label_kwargs : dict
            Any additional keyword arguments passed to matplotlib's Text class for labels and titles
            Used to change font size, color, font family, etc.

        """

        # TODO: Right now you can not have a different number of colorbars on different rows. Row 1 colorbar layout
        #    must be the same as row N.

        # Dark mode
        if dark_mode:
            plt.style.use('dark_background')

        # Initialize various flags
        self._has_suptitle = False

        # nrows and ncols must be integers of at least 1.
        nrows = int(nrows)
        ncols = int(ncols)
        if nrows < 1:
            raise ValueError
        if ncols < 1:
            raise ValueError
        self._nrows = nrows
        self._ncols = ncols

        # Make list of references to the axes.
        self.axis_keys = list()
        for row_i in range(self.nrows):
            for col_i in range(self.ncols):
                self.axis_keys.append((row_i, col_i))

        real_nrows = self.nrows
        self._nsubplots = self.nrows * self.ncols

        # Clean up the colorbar input
        make_colorbars = self._convert_input_to_dict(make_colorbars, input_type='bool')

        # Create width ratio array
        if make_colorbars:
            num_colorbars = sum(list(make_colorbars.values()))

            # Calculate the width of the subplots and colorbars based on the number of colorbars
            individual_width = 1. / (self.ncols + num_colorbars * colorbar_ratio)
            colorbar_width = colorbar_ratio * individual_width

            # Increase the number of "real" axes to include the colorbar axes.
            real_ncols = self.ncols + num_colorbars
        else:
            individual_width = 1. / self.ncols
            colorbar_width = None
            real_ncols = self.ncols

        width_ratios = list()
        for col_i in range(self.ncols):
            # Add the width for each subplot
            width_ratios.append(individual_width)

            # Add the widths for each colorbar after the main axis (if requested by the user)
            index = (0, col_i)
            if index in make_colorbars:
                if make_colorbars[index]:
                    width_ratios.append(colorbar_width)

        # Build gridspec
        figure_size = (figure_scale * aspect_ratio * 2. * self.ncols, figure_scale * 2. * self.nrows)
        self._figure = plt.figure(figsize=figure_size)
        self._gridspec = GridSpec(
            real_nrows, real_ncols, figure=self.figure, width_ratios=width_ratios,
            hspace=hspace, wspace=wspace
            )

        # Determine if any projections are used.
        if projections is None:
            # No projections provided. Assume rectilinear for all plots.
            projections = [None] * self.nsubplots
        else:
            if not CARTOPY_INSTALLED:
                raise ImportError('Projection provided but Cartopy not installed.')

            if type(projections) is str:
                # One projection provided. Assume it is the same for all subplots
                if not projections.lower() in KNOWN_PROJECTIONS:
                    raise KeyError(f'Unknown projection provided for GridPlot: {projections}.')
                projections = [projections] * self.nsubplots
            elif type(projections) is list:
                # Assume that one projection provided for each subplot
                if len(projections) != self.nsubplots:
                    raise ValueError('Unexpected number of projections provided given the number of subplots.')
                try:
                    ['rectilinear' if projection_name is None else
                     KNOWN_PROJECTIONS[projection_name.lower()] for projection_name in projections]
                except KeyError:
                    raise KeyError('One or more unknown projections.')
            else:
                raise AttributeError('Unexpected attributed provided for GridPlot projections.')

            projections = ['rectilinear' if projection_name is None else
                           KNOWN_PROJECTIONS[projection_name.lower()] for projection_name in projections]

        # Build and add axes
        self.axes_by_rowcol = dict()
        self.cb_axes_by_rowcol = dict()
        self.axes_references = dict()
        real_axes = list()
        cb_axes = list()
        real_col_i = 0
        subplot_num = 0
        for col_i in range(self.ncols):

            for row_i in range(self.nrows):
                # Add subplot from gridspec to figure

                # Determine if there are any projections for this subplot
                if projections is not None:
                    projection = projections[subplot_num]
                    if projection is not None:
                        if type(projection) is not str:
                            projection = projection()
                        ax = self.figure.add_subplot(
                            self.gridspec[row_i, real_col_i],
                            projection=projection
                            )
                    else:
                        ax = self.figure.add_subplot(self.gridspec[row_i, real_col_i])
                else:
                    ax = self.figure.add_subplot(self.gridspec[row_i, real_col_i])
                real_axes.append(ax)
                self.axes_by_rowcol[row_i, col_i] = ax
                subplot_num += 1

                # Make colorbars
                index = (0, col_i)
                if index in make_colorbars:
                    if make_colorbars[index]:
                        # Move the real column index forward 1 since we are pulling in the colorbar axis.
                        real_col_i += 1
                        cb_ax = self.figure.add_subplot(self.gridspec[row_i, real_col_i])
                        cb_axes.append(cb_ax)
                        self.cb_axes_by_rowcol[row_i, col_i] = cb_ax

            real_col_i += 1

        self.axes = tuple(real_axes)
        self.cb_axes = tuple(cb_axes)

        # Set subplot scales and labels if provided.
        if label_kwargs is None:
            label_kwargs = dict()
        self.set(
            suptitle=suptitle, titles=titles,
            xlabels=xlabels, xscales=xscales,
            ylabels=ylabels, yscales=yscales,
            reset_values=False, **label_kwargs
            )

        # Allow axes to reference the GridPlot class instance
        for ax in self.axes:
            ax.tpy_gp = self
        for cb_ax in self.cb_axes:
            cb_ax.tpy_gp = self

        # Setup empty containers
        self.colorbars = list()
        self.colorbar_by_rowcol = dict()

    def set_suptitle(self, suptitle: str, top_buffer: float = 0.85):

        self.figure.suptitle(suptitle)

        # Using the tight_layout option will cause the super title to overlap the data. We can adjust the plot to avoid
        #    this, but we only want to do it once. So if a super title was added before and this method is being
        #    called a second time then we don't want to make the adjustment.
        if not self._has_suptitle:
            self.fig.subplots_adjust(top=top_buffer)
            self.fig.subplots_adjust(bottom=0.15)

        # Now we can say that a super title has been added.
        self._has_suptitle = True

        return self

    def _convert_input_to_dict(
        self, inputs: AxisReferencedInput,
        row_priority: bool = False,
        col_priority: bool = False,
        bottom_row_col_priority: bool = False,
        input_type: str = 'str'
        ) -> AxisReferencedOutput:
        """ The user may offer a variety of input structures for various GridPlot methods. This function cleans the
        input and returns a dictionary with axis indices as keys.

        Parameters
        ----------
        inputs : AxisReferencedInput
            User input as either a string, list of strings, tuple of strings, or dictionary of axis
             indices and strings.
        row_priority : bool = False
            If True, inputs will be put across the first column, one for each row.
        col_priority : bool = False
            If True, inputs will be put across the first row, one for each column.
        bottom_row_col_priority : bool = False
            If True, inputs will be put across the bottom row, one for each column.
        input_type : str = 'str'
            Various types of inputs may be passed to this method such as axis labels (strings), axis limits (tuples
                of floats), or booleans. This flag marks what the input, and expected output, types should be.

        Returns
        -------
        cleaned_input : AxisReferencedOutput
            Clean input stored with axis indices as keys.

        """

        # Figure out what kind of input the user is passing in.
        output_type = {'str': str, 'limit': tuple, 'bool': bool}[input_type.lower()]

        cleaned_input = None
        if type(inputs) is str and output_type is str:

            # Check if the input is a string in which case all axes will get the same treatment.
            cleaned_input = {axis_key: inputs for axis_key in self.axis_keys}
        elif type(inputs) is bool and output_type is bool:

            # Check if the input is a bool in which case all axes will get the same treatment.
            cleaned_input = {axis_key: inputs for axis_key in self.axis_keys}
        elif type(inputs) in [dict]:

            # If input is a dict then assume it is already in the correct format, but perform some sanity checks.
            key_0 = list(inputs.keys())[0]
            assert type(key_0) == tuple
            assert type(key_0[0]) == int
            assert type(key_0[1]) == int
            assert type(inputs[key_0]) == output_type
            cleaned_input = inputs
        elif type(inputs) in [list, tuple]:

            if output_type is tuple and type(inputs[0]) is float and len(inputs) == 2:

                # This is a limit input that can be applied to all axes equally.
                cleaned_input = {axis_key: inputs for axis_key in self.axis_keys}
            else:

                # If the input is a list or tuple then it will depend upon the length and assumptions on what to do.
                if len(inputs) == self.nsubplots:

                    # Assume input is in order of subplots
                    cleaned_input = {axis_key: input_str for axis_key, input_str in zip(self.axis_keys, inputs)}
                elif (len(inputs) == self.nrows) and row_priority:

                    # If row priority and the number of inputs equal the number of rows then put the input on the first
                    #   axis of each row.
                    cleaned_input = dict()
                    for row_i in range(self.nrows):
                        cleaned_input[row_i, 0] = inputs[row_i]

                elif (len(inputs) == self.ncols) and col_priority:

                    # If column priority and the number of inputs equal the number of columns then put the input on
                    #    the first axis of each column.
                    cleaned_input = dict()
                    for col_i in range(self.ncols):
                        cleaned_input[0, col_i] = inputs[col_i]

                elif (len(inputs) == self.ncols) and bottom_row_col_priority:

                    # If column priority and the number of inputs equal the number of columns then put the input on
                    #    the first axis of each column.
                    cleaned_input = dict()
                    for col_i in range(self.ncols):
                        cleaned_input[self.nrows - 1, col_i] = inputs[col_i]

        else:
            raise TypeError(f'Unexpected type encountered in GridPlot._convert_input_to_dict: {type(inputs)}.')

        if cleaned_input is None:
            raise AttributeError(f'GridPlot._convert_input_to_dict can not convert input to axis dict: {inputs}.')

        return cleaned_input

    def set_references(self, references: AxisReferencedInputStr, reset_values: bool = False):
        """ Set new references for selected axes.

        Parameters
        ----------
        references : AxisReferencedInput
            New references for selected axes.
        reset_values : bool = False
            If True, then all the axis references will be reset before adding the new ones provided.

        Returns
        -------
        self

        """

        # Clean input
        if references is None:
            return self

        references = self._convert_input_to_dict(references, col_priority=True)

        if reset_values:
            # Clear all user defined references
            self.axes_references = dict()

        for axis_index, new_reference in references.items():

            # Check that the reference is okay. It can not be something already in the GridPlot instance's attributes
            if new_reference in dir(self):
                raise AttributeError(f'New axis reference can not be named after an attribute in GridPlot instance.')

            # Get desired axis
            ax = self.axes_by_rowcol[axis_index]

            # Add new reference to desired axis to axes reference dictionary
            self.axes_references[new_reference] = ax

            # Add new reference to desired axis to instance attributes
            setattr(self, new_reference, ax)

        return self

    def set_titles(self, titles: AxisReferencedInputStr, reset_values: bool = False, **text_kwargs):
        """ Set axis titles for selected subplots.

        Parameters
        ----------
        titles : AxesReferencedInput
            New subplot titles for selected subplots.
        reset_values : bool = False
            If True, then all axes titles will be cleared before they are updated.
        text_kwargs : dict = None
            Optional text or font matplotlib parameters.

        Returns
        -------
        self

        """

        # Clean input
        if titles is None:
            return self

        titles = self._convert_input_to_dict(titles, col_priority=True)

        for axis_index, ax in self.axes_by_rowcol.items():

            if reset_values:
                # Clear titles from all axes
                ax.set_title('')

            if axis_index in titles:
                # Set title for the user requested axes.
                ax.set_title(titles[axis_index], **text_kwargs)

        return self

    def set_xlabels(self, xlabels: AxisReferencedInputStr, reset_values: bool = False, **text_kwargs):
        """ Set x-axis labels for selected subplots.

        Parameters
        ----------
        xlabels : AxesReferencedInput
            New x-axis labels for selected subplots.
        reset_values : bool = False
            If True, then all axes titles will be cleared before they are updated.
        text_kwargs : dict = None
            Optional text or font matplotlib parameters.

        Returns
        -------
        self

        """

        # Clean input
        if xlabels is None:
            return self

        xlabels = self._convert_input_to_dict(xlabels, bottom_row_col_priority=True)

        for axis_index, ax in self.axes_by_rowcol.items():

            if reset_values:
                # Clear x-axis labels from all axes
                ax.set_xlabel('')

            if axis_index in xlabels:
                # Set title for the user requested axes.
                ax.set_xlabel(xlabels[axis_index], **text_kwargs)

        return self

    def set_xscales(self, xscales: AxisReferencedInputStr, reset_values: bool = False):
        """ Set x-axis scales for selected subplots.

        Parameters
        ----------
        xscales : AxesReferencedInput
            New x-axis scales for selected subplots.
        reset_values : bool = False
            If True, then all axes titles will be cleared before they are updated.

        Returns
        -------
        self

        """

        # Clean input
        if xscales is None:
            return self

        xscales = self._convert_input_to_dict(xscales, bottom_row_col_priority=True)

        for axis_index, ax in self.axes_by_rowcol.items():

            if reset_values:
                # Set x scale to linear for all axes
                ax.set_xscale('linear')

            if axis_index in xscales:
                # Set title for the user requested axes.
                ax.set_xscale(xscales[axis_index])

        return self

    def set_ylabels(self, ylabels: AxisReferencedInputStr, reset_values: bool = False, **text_kwargs):
        """ Set y-axis labels for selected subplots.

        Parameters
        ----------
        ylabels : AxesReferencedInput
            New y-axis labels for selected subplots.
        reset_values : bool = False
            If True, then all axes titles will be cleared before they are updated.
        text_kwargs : dict = None
            Optional text or font matplotlib parameters.

        Returns
        -------
        self

        """

        # Clean input
        if ylabels is None:
            return self

        ylabels = self._convert_input_to_dict(ylabels, row_priority=True)

        for axis_index, ax in self.axes_by_rowcol.items():

            if reset_values:
                # Clear y-axis labels from all axes
                ax.set_ylabel('')

            if axis_index in ylabels:
                # Set title for the user requested axes.
                ax.set_ylabel(ylabels[axis_index], **text_kwargs)

        return self

    def set_yscales(self, yscales: AxisReferencedInputStr, reset_values: bool = False):
        """ Set y-axis scales for selected subplots.

        Parameters
        ----------
        yscales : AxesReferencedInput
            New y-axis scales for selected subplots.
        reset_values : bool = False
            If True, then all axes titles will be cleared before they are updated.

        Returns
        -------
        self

        """

        # Clean input
        if yscales is None:
            return self

        yscales = self._convert_input_to_dict(yscales, row_priority=True)

        for axis_index, ax in self.axes_by_rowcol.items():

            if reset_values:
                # Set y scale to linear for all axes
                ax.set_yscale('linear')

            if axis_index in yscales:
                # Set title for the user requested axes.
                ax.set_yscale(yscales[axis_index])

        return self

    def set(
        self, suptitle: str = None, titles: AxisReferencedInputStr = None,
        xlabels: AxisReferencedInputStr = None, xscales: AxisReferencedInputStr = None,
        ylabels: AxisReferencedInputStr = None, yscales: AxisReferencedInputStr = None,
        reset_values: bool = False, **text_kwargs
        ):
        """ Set various labels or scales for selected subplots

        Parameters
        ----------
        suptitle : str = None
            New super title for the figure.
        titles : AxesReferencedInput = None
            New subplot titles for selected subplots.
        xlabels : AxesReferencedInput = None
            New x-axis labels for selected subplots.
        xscales : AxesReferencedInput = None
            New x-axis scales for selected subplots.
        ylabels : AxesReferencedInput = None
            New y-axis labels for selected subplots.
        yscales : AxesReferencedInput = None
            New y-axis scales for selected subplots.
        reset_values : bool = False
            If True, all values will be reset.
        text_kwargs : dict = None
            Optional text or font matplotlib parameters.

        Returns
        -------
        self

        """

        if suptitle is not None:
            self.set_suptitle(suptitle)

        if titles is not None:
            self.set_titles(titles, reset_values=reset_values, **text_kwargs)

        if xlabels is not None:
            self.set_xlabels(xlabels, reset_values=reset_values, **text_kwargs)

        if xscales is not None:
            self.set_xscales(xscales, reset_values=reset_values)

        if ylabels is not None:
            self.set_ylabels(ylabels, reset_values=reset_values, **text_kwargs)

        if yscales is not None:
            self.set_yscales(yscales, reset_values=reset_values)

    def set_colorbar_data(self, cb_data, cb_label: str = None):

        # Clear any old data (if set_colorbar_data was called before)
        if self.colorbars != list():
            self.colorbars = list()

        if self.colorbar_by_rowcol != dict():
            self.colorbar_by_rowcol = dict()

        # Make a colorbar for each colorbar axis
        for (row_i, col_i), colorbar_ax in self.cb_axes_by_rowcol.items():
            colorbar = plt.colorbar(cb_data, cax=colorbar_ax)

            if cb_label is not None:
                colorbar.set_label(cb_label)

            self.colorbars.append(colorbar)
            self.colorbar_by_rowcol[row_i, col_i] = colorbar

    def savefig(self, filename: str, save_png: bool = False, png_dpi: int = 300):

        #  Find a unique figure path name
        if '.pdf' == filename[-4:]:
            filename = filename[:-4]
        filename = unique_path(filename + '.pdf', is_dir=False)
        filepath = filename[:-4]

        # Save figure pdf (scalable graphics)
        self.figure.savefig(filepath + '.pdf')

        # Save figure png
        if save_png:
            self.figure.savefig(filepath + '.png', dpi=png_dpi)

    @staticmethod
    def show():
        plt.show()

    # # Private properties
    @property
    def ncols(self) -> int:
        """ Return the number of columns in GridPlot. """
        return self._ncols

    @ncols.setter
    def ncols(self, value):
        raise AttributeError('Property `ncols` can not be changed after GridPlot instance has been created.')

    @property
    def nrows(self) -> int:
        """ Return the number of rows in GridPlot. """
        return self._nrows

    @nrows.setter
    def nrows(self, value):
        raise AttributeError('Property `nrows` can not be changed after GridPlot instance has been created.')

    @property
    def nsubplots(self) -> int:
        """ Return the number of subplots in GridPlot. """
        return self._nsubplots

    @nsubplots.setter
    def nsubplots(self, value):
        raise AttributeError('Property `nsubplots` can not be changed after GridPlot instance has been created.')

    @property
    def gridspec(self) -> GridSpec:
        """ Return the GridSpec instance used in GridPlot. """
        return self._gridspec

    @gridspec.setter
    def gridspec(self, value):
        raise AttributeError('Property `gridspec` can not be changed after GridPlot instance has been created.')

    @property
    def figure(self) -> plt.Figure:
        """ Return the matplotlib.Figure instance used in GridPlot. """
        return self._figure

    @figure.setter
    def figure(self, value):
        raise AttributeError('Property `figure` can not be changed after GridPlot instance has been created.')

    # # Aliased properties
    @property
    def gs(self) -> GridSpec:
        """ Return the GridSpec instance used in GridPlot. """
        return self.gridspec

    @gs.setter
    def gs(self, value):
        self.gridspec = value

    @property
    def fig(self) -> plt.Figure:
        """ Return the matplotlib.Figure instance used in GridPlot. """
        return self.figure

    @fig.setter
    def fig(self, value):
        self.figure = value

    # # Dunder methods
    def items(self):
        return self.axes_by_rowcol.items()

    def __getitem__(self, key: Union[str, Tuple[int, int]]) -> plt.Axes:

        # First check if the key is in the user defined reference dictionary.
        if key in self.axes_references:
            return self.axes_references[key]

        # Next try the index dictionary
        if key in self.axes_by_rowcol:
            return self.axes_by_rowcol[key]

        # No luck :(
        raise KeyError

    def __setitem__(self, key, value):

        self.axes_references[key] = value

    def __iter__(self):
        return iter(self.axes)
