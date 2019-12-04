from .. import version
from .dictionary_utils import nested_get
from .logging import TidalLogger
from .progress import progress_bar


class TidalPyClass:
    """ All functional classes used in TidalPy inherit from this class
    """

    tidalpy_version = version

    def __init__(self):
        self.pyname = f'{self.__class__.__name__}'