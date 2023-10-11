def test_version():
    # Test Load
    import TidalPy
    TidalPy.test_mode()

    from TidalPy import version

    assert version is not None
