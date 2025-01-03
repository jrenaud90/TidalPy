import operator
from typing import Callable, TYPE_CHECKING, Tuple, Union

import TidalPy
from TidalPy.exceptions import (AttributeNotSetError, ConfigPropertyChangeError, InitiatedPropertyChangeError,
                                MissingArgumentError, OuterscopePropertySetError, ParameterMissingError,
                                UnknownModelError)

from TidalPy.utilities.dictionary_utils import nested_get, nested_place
from TidalPy.utilities.classes.config.config import ConfigHolder

from TidalPy.logger import get_logger
log = get_logger("TidalPy")

if TYPE_CHECKING:
    from TidalPy.structures.layers import PhysicalLayerType
    from TidalPy.structures.world_types import LayeredWorldType
    from TidalPy.utilities.types import NoneType


class ModelHolder(ConfigHolder):
    """ ModelHolder
    This class serves as a parent for most physics models (e.g., Radiogenics, Rheology, etc.)

    Provides basic functionality to load in default model inputs and run calculations using those inputs.
    """

    known_models = None
    known_model_const_args = None
    known_model_live_args = None

    def __init__(self, model_name: str = None, replacement_config: dict = None, initialize: bool = True):
        """ Constructor for ModelHolder class

        Parameters
        ----------
        model_name : str = None
            The user-provided model name.
        replacement_config : dict = None
            Optional dictionary that will replace TidalPy's default configurations for a particular model.
        initialize : bool = True
            Determines if initial reinit should be performed on the model (loading in data from its `self.config`).
        """

        # Setup parent class
        super().__init__(replacement_config=replacement_config)

        # Pull out model information, check if it is a valid model name, then store it
        if model_name is None and 'model' in self.config:
            model_name = self.config['model']

        # Configuration properties
        self._model = model_name
        self._func = None
        # Some functions can support arrays and floats interchangeably. Others require separate functions.
        #    These functions should be stored separately.
        self._func_array = None
        # If a separately defined function was defined then set this flag to True.
        self._func_array_defined = False
        # Functions may have separate inputs. These are stored as either:
        #    `inputs`: tuple of constants passed to the self.func or self.func_array
        #    `live_inputs`: tuple of dynamic inputs that may change after TidalPy is loaded
        #       (e.g., the viscosity of a layer)
        self.get_live_args = None
        self._live_inputs = None
        self._inputs = None
        self._constant_arg_names = None
        self._live_arg_names = None
        # Calculation properties and methods
        self._debug_mode_on = False
        self._calc_to_use = None  # type: Union[None, Callable]

        if initialize:
            self.reinit(initial_init=True)

    def reinit(self, initial_init: bool = False):
        """ Model will look at the user-provided configurations and pull out model information including constants

        Parameters
        ----------
        initial_init : bool = False
            Must be set to `True` if this is the first time this method has been called (additional steps may be
                preformed during the first reinit call).
        """

        if initial_init:
            log.debug(f'Initializing {self}.')
        else:
            log.debug(f'Reinit called for {self}.')

        # Build constant and live inputs, and reinit self.func
        self.build_inputs()

        # Switch between calculate and calculate debug.
        #    Generally speaking, _calculate_debug is a much slower function that includes additional sanity checks.
        if TidalPy.extensive_checks:
            if '_calculate_debug' in self.__dict__:
                self._calc_to_use = getattr(self, '_calculate_debug')
                self._debug_mode_on = True
            else:
                log.debug(
                    f"TidalPy's debug mode is on, but no debug calculation method found for {self}. "
                    f"Using Regular."
                    )
                self._calc_to_use = getattr(self, '_calculate')
        else:
            self._calc_to_use = getattr(self, '_calculate')

        if self._calc_to_use is None:
            raise AttributeNotSetError

        # TODO: The below breaks when calculate is a method. Not sure if there is a solution:
        #    Error: AttributeError: attribute '__doc__' of 'method' objects is not writable
        # # Give calculate the same doc string as whatever is store in _calc (debug or regular)
        # if self._calc.__doc__ not in [None, '']:
        #     self.calculate.__doc__ = self._calc.__doc__

    def calculate(self, *args, **kwargs):
        """ Main calculation point for the model """

        if self._calc_to_use is None:
            raise AttributeNotSetError

        # Some models have inputs that need to be updated at each call
        if self.get_live_args is not None:
            try:
                self._live_inputs = self.get_live_args()
            except AttributeError:
                raise AttributeNotSetError('One or more live arguments are not set or their references were not found.')

        return self._calc_to_use(*args, **kwargs)

    def build_inputs(self):
        """ Build a list of the live and constant inputs for the model's calculate function.
        """

        # Try to find the provided model and load in its main function
        test_names = [self.model, self.model.lower(), self.model.title()]
        model_found = False
        model_name = None
        for test_name in test_names:
            if test_name in self.known_models:
                self._func = self.known_models[test_name]
                model_found = True
                model_name = test_name
                break

        if not model_found:
            raise UnknownModelError(f'Unknown model: "{self.model}" for {self}')
        elif model_name != self.model:
            log.warning(
                f'Model "{self.model}" for {self} does not match TidalPy capitalization. '
                f'Changing to {model_name}.'
                )
            self._model = model_name

        # Load the model's functions
        if self.model + '_array' in self.known_models:
            self._func_array = self.known_models[self.model + '_array']
            self._func_array_defined = True
        else:
            # No separate array function found. Try to use the regular float version
            self._func_array = self.func
            self._func_array_defined = False

        # Determine what (if any) constant and live arguments exist for this model's functions
        self._constant_arg_names = self.known_model_const_args[self.model]
        self._live_arg_names = self.known_model_live_args[self.model]

        # Pull out constant arguments
        self._inputs = self.build_args(self._constant_arg_names, parameter_dict=self.config, is_live_args=False)

        # Build live argument functions and calls
        self._live_inputs = None
        if len(self._live_arg_names) == 0:
            self.get_live_args = None
        else:
            live_funcs = self.build_args(self._live_arg_names, is_live_args=True)
            # The live_funcs is a list of python operator.attrgetter functions that will pull a property from whatever
            #    is called to them. For example the func for 'time' when you do `func(self)` will pull `self.time`
            self.get_live_args = lambda: tuple([live_func(self) for live_func in live_funcs])

    # # Configuration properties
    @property
    def model(self) -> str:
        """ Model name """
        return self._model

    @model.setter
    def model(self, value):
        raise ConfigPropertyChangeError

    @property
    def func(self) -> Callable:
        """ Callable function based on the user-provided model.

        See Also
        --------
        ModelHolder.func_array
        """

        return self._func

    @func.setter
    def func(self, value):
        raise ConfigPropertyChangeError

    @property
    def func_array(self) -> Callable:
        """ Callable function based on the user-provided model.

         Notes
         -----
         .. If `func_array_defined` is True this this will be a different callable from `self.func` designed to
            specifically work with numpy arrays. Otherwise, the `func` and `func_array` are identical.

        See Also
        --------
        ModelHolder.func
        """

        return self._func_array

    @func_array.setter
    def func_array(self, value):
        raise ConfigPropertyChangeError

    @property
    def func_array_defined(self) -> bool:
        """ Flag for if the `func_array` property is set or not """
        return self._func_array_defined

    @func_array_defined.setter
    def func_array_defined(self, value):
        raise ConfigPropertyChangeError

    @property
    def inputs(self) -> Union['NoneType', Tuple[float, ...]]:
        """ Some models may require additional constants to be passed to the `self.func` or `self.func_array`.
            These are stored in this tuple if applicable.
        """

        return self._inputs

    @inputs.setter
    def inputs(self, value):
        raise ConfigPropertyChangeError

    @property
    def live_inputs(self) -> Union['NoneType', Tuple[float, ...]]:
        """ Similar to `self.inputs` but these are dynamic parameters that can change after initialization
            (e.g., the viscosity of a layer).
        """

        return self._live_inputs

    @live_inputs.setter
    def live_inputs(self, value):
        raise ConfigPropertyChangeError

    # Calculation properties
    @property
    def debug_mode_on(self) -> bool:
        """ Flag for if the model's debug function is being used """
        return self._debug_mode_on

    @debug_mode_on.setter
    def debug_mode_on(self, value):
        raise ConfigPropertyChangeError

    # # Dunder properties
    def __str__(self):
        name_str = f'{self.__class__.__name__}'
        if '_model' in self.__dict__:
            if self.model is not None:
                name_str += ' ({self.model})'
        return name_str

    @staticmethod
    def build_args(arg_names: Tuple[str, ...], parameter_dict: dict = None, is_live_args: bool = False):
        """ Build an input tuple based on the required, constant, arguments and a parameter dictionary

        Parameters
        ----------
        arg_names : Tuple[str, ...]
            List of required constant argument names needed for this model's function.
        parameter_dict : dict
            Dictionary of parameters.
        is_live_args : bool
            Flag if this call is looking for live args, rather than constant args.

        Returns
        -------
        args : Tuple[Any, ...]
            List of default argument parameters
        """

        args = list()

        if not is_live_args:
            # Constant Arguments
            if parameter_dict is None:
                raise MissingArgumentError(
                    f'Parameter configuration dictionary is required to build constant '
                    f'model arguments.'
                    )

            for arg_name in arg_names:
                if arg_name not in parameter_dict:
                    raise ParameterMissingError(f'Parameter: {arg_name} is missing from configuration dictionary.')
                args.append(parameter_dict[arg_name])
        else:
            # Live Arguments
            for live_arg_signature in arg_names:

                # Create a attrgetter for the live argument.
                if 'self.' in live_arg_signature:
                    live_arg_signature = live_arg_signature.split('self.')[1]

                getter_func = operator.attrgetter(live_arg_signature)
                # The getter_func is a python operator.attrgetter function that will pull a property from whatever
                #    is called to them. For example the func for 'time' when you do `func(self)` will pull `self.time`
                args.append(getter_func)

        return tuple(args)


class LayerModelHolder(ModelHolder):
    """ LayerModelHolder
    Parent class for physics models that are stored within a world's layer

    Provides basic functionality to load in default model inputs and run calculations using those inputs and the layer's
        current state properties (e.g., temperature).

    See Also
    --------
    TidalPy.cooling
    TidalPy.radiogenics
    TidalPy.rheology
    """

    model_config_key = None

    def __init__(
        self, layer: 'PhysicalLayerType', model_name: str = None, store_config_in_layer: bool = True,
        initialize: bool = True
        ):
        """ Constructor for LayerModelHolder class

        Parameters
        ----------
        layer : PhysicalLayerType
            The layer instance which the model should perform calculations on.
        model_name : str = None
            The user-provided model name.
        store_config_in_layer: bool = True
            Flag that determines if the final model's configuration dictionary should be copied into the
            `layer.config` dictionary.
        initialize : bool = True
            Determines if initial reinit should be performed on the model (loading in data from its `self.config`).
        """

        # Initialized properties
        # Store layer and world information
        self._layer = layer
        self._world = layer.world

        # Get default configurations for the model
        self.default_config = TidalPy.config['layers'][layer.type][self.model_config_key]

        # Record if model config should be stored back into layer's config
        self._store_config_in_layer = store_config_in_layer

        # Find the model's new configurations which the user is expected to have stored in the layer's configuration
        #    dictionary.

        # The layer's type is used to pull out default parameter information
        self.default_config_key = self.layer_type

        config = None
        try:
            config = nested_get(self.layer.config, self.model_config_key, raiseon_nolocate=True)
        except KeyError:
            log.debug(f"User provided no model information for {self}; using defaults instead.")

        if config is None and self.default_config is None:
            raise ParameterMissingError(f"Config not provided for {self}; and no defaults are set.")

        # Setup ModelHolder and ConfigHolder methods. Using the layer's config file as the replacement config.
        super().__init__(model_name=model_name, replacement_config=config, initialize=initialize)

        if self.store_config_in_layer:
            # Once the configuration file is constructed (with defaults and any user-provided replacements) then
            #    store the new config in the layer's config, overwriting any previous parameters.
            nested_place(self.config, self.layer.config, self.model_config_key, make_copy=False, retain_old_value=True)

    # # Initialized properties
    @property
    def layer(self) -> 'PhysicalLayerType':
        """ Layer instance which the model performs calculations on """
        return self._layer

    @layer.setter
    def layer(self, value):
        raise InitiatedPropertyChangeError

    @property
    def world(self) -> 'LayeredWorldType':
        """ Layered world where the `self.layer` is stored """
        return self._world

    @world.setter
    def world(self, value):
        raise InitiatedPropertyChangeError

    @property
    def store_config_in_layer(self) -> bool:
        """ Flag to store model's configuration dictionary in the `layer.config` dictionary """
        return self._store_config_in_layer

    @store_config_in_layer.setter
    def store_config_in_layer(self, value):
        raise InitiatedPropertyChangeError

    # # Outer-scope Properties
    @property
    def layer_type(self):
        return self.layer.type

    @layer_type.setter
    def layer_type(self, value):
        raise OuterscopePropertySetError

    # # Dunder properties
    def __str__(self):
        name_str = f'{self.__class__.__name__}'
        if '_model' in self.__dict__:
            if self.model is not None:
                name_str += f' ({self.model}'
        if '_layer' in self.__dict__:
            if self.layer is not None:
                if '(' in name_str:
                    name_str += f' in {self.layer.name})'
                else:
                    name_str += f' (Layer: {self.layer.name})'

        if '(' in name_str and ')' not in name_str:
            name_str += ')'

        return name_str
