burnman_installed = True
try:
    from burnman.classes.material import Material, material_property
    from burnman.classes.mineral import Mineral
    from burnman.tools.chemistry import dictionarize_formula, formula_mass
except ImportError:
    burnman_installed = False
    # Build fake class so type checking passes.
    class burnman:
        def __init__(self):
            Planet = None
            Layer = None

    class Mineral:
        pass

    class Material:
        pass

    dictionarize_formula = lambda x: None
    formula_mass = lambda x: None
    material_property = lambda x: None


from .ice import *
from .pyrite import Pyrite
