from ...version import version


class TidalPyClass:
    """ All functional classes used in TidalPy inherit from this class
    """

    tidalpy_version = version

    def __init__(self):
        self.pyname = f'{self.__class__.__name__}'