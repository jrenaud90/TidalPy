import copy
from typing import Tuple, TYPE_CHECKING
import operator

from .classes import ConfigHolder
from .. import debug_mode
from ..exceptions import (ImplementedBySubclassError, MissingArgumentError, ParameterMissingError, TidalPyException,
                          UnknownModelError, AttributeNotSetError)
from ..initialize import log

if TYPE_CHECKING:
    from ..structures.layers import ThermalLayer


class ModelHolder(ConfigHolder):

    """ Classes which are used in OOP calculation scheme
    """

    known_models = None
    known_model_const_args = None
    known_model_live_args = None

    def __init__(self, model_name: str = None, replacement_config: dict = None, call_reinit: bool = True):

        # Pull out model information, check if it is a valid model name, then store it
        if model_name is None and 'model' in self.config:
            model_name = self.config['model']

        # Store model name information
        self.model = model_name
        self.pyname += '_' + self.model

        super().__init__(replacement_config=replacement_config, call_reinit=call_reinit)

        # Attempt to find the model's function information
        self.func = None
        self._constant_arg_names = None
        self._live_arg_names = None
        try:
            self.func = self.known_models[self.model]
            self._constant_arg_names = self.known_model_const_args[self.model]
            self._live_arg_names = self.known_model_live_args[self.model]
        except KeyError:
            raise UnknownModelError(f'Unknown model: {self.model} for {self.__class__.__name__}')

        # Pull out constant arguments
        self.inputs = self.build_args(self._constant_arg_names, parameter_dict=self.config, is_live_args=False)

        # Build live argument functions and calls
        self.live_inputs = None
        if len(self._live_arg_names) == 0:
            self.get_live_args = None
        else:
            live_funcs = self.build_args(self._live_arg_names, is_live_args=True)
            self.get_live_args = lambda: tuple([live_func(self) for live_func in live_funcs])

        # Switch between calculate and calculate debug. Generally, _calculate_debug is a much slower function that
        #    includes additional checks
        if debug_mode:
            if '_calculate_debug' in self.__dict__:
                self._calc = self._calculate_debug
            else:
                log(f'Debug mode is on, but it appears that no debug calculation method has been implemented '
                    f'for {self.__class__.__name__}. Using regular calculate method.')
                self._calc = self._calculate
        else:
            self._calc = self._calculate

    def calculate(self, *args, **kwargs):

        # Some models have inputs that need to be updated at each call
        if self.get_live_args is not None:
            try:
                self.live_inputs = self.get_live_args()
            except AttributeError:
                raise AttributeNotSetError('One or more live arguments are not set.')

        return self._calc(*args, **kwargs)

    def _calculate(self, *args, **kwargs):
        raise ImplementedBySubclassError

    # def _calculate_debug(self, *args, other_inputs: tuple = tuple(), **kwargs):
    #     raise ImplementedBySubclassError

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
                raise MissingArgumentError('Parameter configuration dictionary is required to build constant args.')

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
                args.append(getter_func)

        return tuple(args)


class LayerModelHolder(ModelHolder):

    """ Models which are stored within a layer and make calls to that layer's attributes and methods.
    """

    model_config_key = None

    def __init__(self, layer: ThermalLayer, model_name: str = None,
                 store_config_in_layer: bool = True, call_reinit: bool = True):

        # Store layer and world information
        self.layer = layer
        self.world = None
        self.layer_type = layer.type
        world_name = 'Unknown'
        if 'world' in layer.__dict__:
            self.world = layer.world
            world_name = self.world.name

        # Update pyname
        self.pyname += '_' + self.layer_type

        # The layer's type is used to pull out default parameter information
        self.default_config_key = self.layer_type

        config = None
        try:
            if type(self.model_config_key) == str:
                config = self.layer.config[self.model_config_key]
            elif type(self.model_config_key) in [list, tuple]:
                config = self.model_config_key[0]
                for subkey in self.model_config_key[1:]:
                    config = config[self.model_config_key]

        except KeyError:
            log(f"User provided no model information for [layer: {self.layer.name} in world: {world_name}]'s "
                f"{self.__class__.__name__}, using defaults instead.", level='debug')

        if config is None and self.default_config is None:
            raise ParameterMissingError(f"Config was not provided for [layer: {self.layer.name} in world: {world_name}]'s "
                                        f"{self.__class__.__name__} and no defaults are set.")

        # Setup ModelHolder and ConfigHolder classes. Using the layer's config file as the replacement config.
        super().__init__(model_name=model_name, replacement_config=config, call_reinit=call_reinit)

        if store_config_in_layer:
            # Once the configuration file is constructed (with defaults and any user-provided replacements) then
            #    store the new config in the layer's config, overwriting any previous parameters.
            if self.model_config_key in self.layer.config:
                # Store the old config under a new key
                self.layer._config[f'OLD_{self.model_config_key}'] = self.layer.config[self.model_config_key]
            self.layer._config[self.model_config_key] = copy.deepcopy(self.config)