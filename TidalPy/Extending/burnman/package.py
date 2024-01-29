import warnings

burnman_installed = True
try:
    import burnman
    from burnman.classes.material import Material, material_property
    from burnman.classes.mineral import Mineral
    from burnman.tools.chemistry import dictionarize_formula, formula_mass
except ImportError:
    warnings.warn("BurnMan installation can not be found. TidalPy's BurnMan extension functions can not be used.")
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

    dictionarize_formula = lambda x: None
    formula_mass = lambda x: None
    material_property = lambda x: None
