import TidalPy
from TidalPy.utilities.dictionary_utils import nested_merge
from TidalPy.logger import get_logger
log = get_logger("TidalPy")

burnman_installed = True
try:
    log.debug('Attempting to import the BurnMan package.')
    import burnman    
except ImportError:
    log.warning("BurnMan installation can not be found. TidalPy's BurnMan extension functions can not be used.")
    burnman_installed = False
    # Build fake class so type checking passes.
    class burnman:
        Planet = None
        Layer = None
        Material = None
    
    class Mineral:
        pass

    class Material:
        pass

    def dictionarize_formula(x):
        return None
    def formula_mass(x):
        return None
    def material_property(x):
        return None

else:
    log.debug(f'BurnMan version {burnman.__version__} was found!')
    from burnman.classes.material import Material as Material
    from burnman.classes.material import material_property as material_property
    from burnman.classes.mineral import Mineral as Mineral
    from burnman.tools.chemistry import dictionarize_formula as dictionarize_formula
    from burnman.tools.chemistry import formula_mass as formula_mass

    log.debug('Appending TidalPy config with BurnMan specific configurations.')
    from .burnman_defaultc import default_burnman_configs

    TidalPy.config = nested_merge(TidalPy.config, default_burnman_configs, make_copies=False)
    log.debug('BurnMan extensions initialized.')
    
    from .burnman_layer import BurnmanLayer as BurnmanLayer
    from .burnman_world import BurnManWorld as BurnManWorld