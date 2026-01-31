import TidalPy

def test_imports():

    tidalpy_includes = TidalPy.get_include()

    assert type(tidalpy_includes) == list
    assert len(tidalpy_includes) > 0


    import CyRK
    cyrk_includes = CyRK.get_include()
    for dir_ in cyrk_includes:
        assert dir_ in tidalpy_includes
