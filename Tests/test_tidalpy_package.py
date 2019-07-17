
def test_import_tidalpy():

    import TidalPy
    # Just do something to ensure TidalPy is loaded into memory and initialized.
    assert type(TidalPy.version) == str
