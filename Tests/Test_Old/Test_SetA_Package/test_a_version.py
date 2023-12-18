def test_version():
    # Test Load
    import TidalPy
    from TidalPy import version

    assert version is not None
