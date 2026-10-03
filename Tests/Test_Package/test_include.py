"""Tests for ``TidalPy.get_include``."""
import os

import CyRK
import pytest

import TidalPy

TIDALPY_DIR = os.path.dirname(os.path.abspath(TidalPy.__file__))


def test_includes_cyrk():
    """The include list holds every CyRK include directory."""
    tidalpy_includes = TidalPy.get_include()
    assert isinstance(tidalpy_includes, list)
    for directory in CyRK.get_include():
        assert directory in tidalpy_includes


@pytest.mark.parametrize(
    "header",
    [
        "constants_.hpp",
        os.path.join("Utilities", "arrays", "interp_.hpp"),
        os.path.join("Utilities", "math", "numerics_.hpp"),
        os.path.join("RadialSolver", "rs_solution_.hpp"),
        os.path.join("Material", "eos", "eos_solution_.hpp"),
    ])
def test_includes_tidalpy_header_directory(header):
    """The directory of each widely included TidalPy header is listed."""
    header_path = os.path.join(TIDALPY_DIR, header)
    assert os.path.isfile(header_path)
    assert os.path.dirname(header_path) in TidalPy.get_include()
