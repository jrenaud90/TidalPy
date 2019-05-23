import burnman

from TidalPy.exceptions import UnknownModelError
from ... import debug_mode

known_materials = None
known_materials_sourceless = None

def setup_material_lists():
    """ Builds a dictionary of known materials

    This is a function, rather than in the body of the script, because it may be a slow process that is only needed
        when the user is building a new planet.
    """

    global known_materials
    global known_materials_sourceless

    # Load in official BurnMan materials
    ignore_list = ['__name__', '__doc__', '__package__', '__loader__', '__spec__', '__path__',
                    '__file__', '__cached__', '__builtins__', 'absolute_import', 'np', 'os', 'jit', 'njit',
                    'logish', 'helpers', 'Mineral', 'SolidSolution']

    # Search burnman.minerals for new material files.
    known_materials_tmp = dict()
    for material_file in burnman.minerals.__dict__:
        if material_file in ignore_list:
            continue
        known_materials_tmp[material_file] = dict()

        # Search a material file for new materials
        for item in burnman.minerals.__dict__[material_file].__dict__:
            if item in ignore_list:
                continue
            potential_obj = burnman.minerals.__dict__[material_file].__dict__[item]
            try:
                # Ignore anything that is not a BurnMan mineral
                if issubclass(potential_obj, burnman.Mineral):
                    known_materials_tmp[material_file][item] = potential_obj
            except TypeError:
                pass

    # Clean-up known materials by lowering all keys
    known_materials = dict()
    known_materials_sourceless = dict()
    for source in known_materials_tmp:
        known_materials[source.lower()] = dict()
        for material in known_materials_tmp[source]:
            known_materials[source.lower()][material.lower()] = known_materials_tmp[source][material]
            if material in known_materials_sourceless:
                # known_materials_sourceless does not know what file a material came from. If 2+ materials with the
                #    same name are found, only the first one will be stored in this dictionary.
                # NOTE: It is always better to provide a material_source filename!!
                pass
            else:
                known_materials_sourceless[material] = (source, known_materials_tmp[source][material])

    # Load in TidalPy custom materials
    # TODO


def find_material(material_name: str, material_source: str = None):
    """ Finds a BurnMan material from a specified source file (if provided)

    :return: returns an un-initialized burnman.material class
    """

    if known_materials is None or known_materials_sourceless is None:
        setup_material_lists()
    
    material_class = None
    if material_source is None:
        try:
            material_source, material_class = known_materials_sourceless[material_name]
        except KeyError:
            UnknownModelError(f'Unknown material: {material_name}. No source filename was provided. Providing a source filename may correct this error.')
    else:
        if material_source not in known_materials:
            raise UnknownModelError(f'Material source filename {material_source} not found. If source file unknown set material_source to None.')
        if material_name not in known_materials[material_source]:
            raise UnknownModelError(f'Material {material_name} not found in source file: {material_source}')
        material_class = known_materials[material_source][material_name]

    return material_class