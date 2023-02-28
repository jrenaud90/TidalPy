from .. import burnman, burnman_installed

if burnman_installed:
    from burnman.classes.material import Material, material_property
    from burnman.classes.mineral import Mineral
    from burnman.tools.chemistry import dictionarize_formula, formula_mass
else:
    burnman_installed = False
    # Build fake classes so type checking passes.

    class Mineral:
        pass

    class Material:
        pass

    dictionarize_formula = lambda x: None
    formula_mass = lambda x: None
    material_property = lambda x: None


from .ice import *
from .pyrite import Pyrite
