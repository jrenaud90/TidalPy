import os

import TidalPy


def test_get_include():
    """`get_include` returns TidalPy's C++ source directories plus every CyRK include directory."""
    import CyRK

    tidalpy_includes = TidalPy.get_include()

    assert type(tidalpy_includes) is list
    assert len(tidalpy_includes) > 0

    for dir_ in CyRK.get_include():
        assert dir_ in tidalpy_includes

    # The TidalPy entries point at directories that ship the C++ headers.
    for dir_ in tidalpy_includes[len(CyRK.get_include()):]:
        assert os.path.isdir(dir_)
        assert any(name.endswith((".hpp", ".cpp")) for name in os.listdir(dir_))
