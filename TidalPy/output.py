from pathlib import Path

import TidalPy
from TidalPy.logger import get_logger

log = get_logger("TidalPy")

def set_output_dir(new_output_dir: str) -> str:
    """Sets new output directory for TidalPy data and logs.

    Parameters
    ----------
    new_output_dir : str
        New output directory. TidalPy will create a new directory if it does not exist.
    
    Returns
    -------
    str
        New output directory.
    """

    assert type(new_output_dir) is str
    log.debug(f'TidalPy output directory changing to {new_output_dir}.')
    TidalPy._output_path = new_output_dir
    return new_output_dir

def create_output_dir() -> str:
    """Creates an output directory for TidalPy data.

    Returns
    -------
    str
        Path to output directory.
    """

    Path(TidalPy._output_path).mkdir(parents=True, exist_ok=True)
    return TidalPy._output_path

