def test_import_tidalpy():
    import TidalPy

    TidalPy.verbose_level = 0
    TidalPy.logging_level = 0
    TidalPy.use_disk = False

    # Just do something to ensure TidalPy is loaded into memory and initialized.
    assert type(TidalPy.version) == str
