def test_version():
    # Test Load
    import TidalPy
    TidalPy.test_mode()

    from TidalPy._version import version

    assert version is not None
