from .layered import LayeredTides

if TYPE_CHECKING:
    from ...structures.world_types import LayeredWorld
    from ...structures.layers import PhysicalLayerType


class MultilayerTides(LayeredTides):
    """ MultilayerTides
    Class used for layered planets (icy or rocky world_types)

    Utilizes shooting or propagation method to solve the viscoelastic-gravitational problem to find a planet's
        Love numbers and radial solutions ($y_{i}$).

    See Also
    --------
    TidalPy.tides.methods.LayeredTides
    """

    model = 'multilayer-base'
    default_config = tide_defaults['layered_love']
    solver_method = 'base'

    def __init__(self, world: 'LayeredWorld', store_config_in_world: bool = True, initialize: bool = True):
        """ Constructor for LayeredLove class

        Parameters
        ----------
        world : TidalWorldType
            The world where tides are being calculated.
        store_config_in_world : bool = True
            Flag that determines if the final model's configuration dictionary should be copied into the
            `world.config` dictionary.
        initialize : bool = True
            If `True`, then an initial call to the tide's reinit method will be made at the end of construction.
        """

        super().__init__(world, store_config_in_world=store_config_in_world, initialize=initialize)

        # State variables
        self._love_k = None
        self._love_h = None
        self._love_l = None
        self._tidal_stress_by = None
        self._tidal_strain_by = None
        self._volumetric_heating = None
        self._total_potential = None
        self._tidal_potential = None
        self._complex_shear = None
        self._tidal_modes =


        heating, volumetric_heating, strains, stresses, total_potential, tidal_potential, \
        complex_shears_avg, tidal_y_avg, love_results, tidal_modes, modes_skipped

    def reinit(self, initial_init: bool = False):
        """ Load configurations into the Tides class and import any config-dependent functions.

        This reinit process is separate from the __init__ method because the Orbit class may need to overload some
            configurations after class initialization.

        Parameters
        ----------
        self
        initial_init : bool = False
            This should be set to True the first time reinit is called.
        """

    def _y_solver(self):
        """ Solve the viscoelastic-gravitational problem to find the y-solutions in the planet """

    # # State properties
    @property
    def tidal_y(self) -> Tuple[float, float, float, float]:
        """ The inputs required to calculate tides - these could change dynamically so they need to be pulled live """
        return self.world.tidal_scale, self.radius, self.world.density_bulk, self.world.gravity_surface

    @tidal_inputs.setter
    def tidal_inputs(self, value):
        raise IncorrectMethodToSetStateProperty