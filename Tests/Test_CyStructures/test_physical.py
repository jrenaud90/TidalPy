import numpy as np
import pytest

from TidalPy.cystructures.physical import PhysicalStructure


def test_physical_structure_radial():
    """Test TidalPy's Physical Structure Cython Class"""
    
    # Test providing nothing (nan'd structure)
    nanlayer = PhysicalStructure("nanlayer")
    assert np.isnan(nanlayer.radius)
    assert np.isnan(nanlayer.radius_upper)
    assert np.isnan(nanlayer.radius_lower)
    assert np.isnan(nanlayer.volume)
    assert np.isnan(nanlayer.surface_area_upper)
    assert np.isnan(nanlayer.surface_area_lower)

    # Test basic physical structure
    layer = PhysicalStructure("This is a layer", radius_upper=10., radius_lower=5.)

    assert layer.name == "This is a layer"
    assert layer.radius == 10.
    assert layer.radius_upper == 10.
    assert layer.radius_lower == 5.
    assert layer.thickness == 5.

    assert layer.volume == (4. / 3.) * np.pi * (10.**3 - 5.**3)
    assert layer.surface_area_upper == 4. * np.pi * 10.**2
    assert layer.surface_area_lower == 4. * np.pi * 5.**2

    # Test providing it a upper radius and thickness instead
    layer = PhysicalStructure("This is a layer", radius_upper=10., thickness=5.)

    assert layer.name == "This is a layer"
    assert layer.radius == 10.
    assert layer.radius_upper == 10.
    assert layer.radius_lower == 5.
    assert layer.thickness == 5.

    assert layer.volume == (4. / 3.) * np.pi * (10.**3 - 5.**3)
    assert layer.surface_area_upper == 4. * np.pi * 10.**2
    assert layer.surface_area_lower == 4. * np.pi * 5.**2

    # Test providing it a lower radius and thickness instead
    layer = PhysicalStructure("This is a layer", radius_lower=10., thickness=5.)

    assert layer.name == "This is a layer"
    assert layer.radius == 15.
    assert layer.radius_upper == 15.
    assert layer.radius_lower == 10.
    assert layer.thickness == 5.

    assert layer.volume == (4. / 3.) * np.pi * (15.**3 - 10.**3)
    assert layer.surface_area_upper == 4. * np.pi * 15.**2
    assert layer.surface_area_lower == 4. * np.pi * 10.**2


def test_physical_structure_radial_error():
    """Test for bad radial states."""
    
    # Provide bad values
    with pytest.raises(ValueError):
        layer = PhysicalStructure("BadLayer", radius_upper=-10.)
    with pytest.raises(ValueError):
        layer = PhysicalStructure("BadLayer", radius_lower=-10.)
    with pytest.raises(ValueError):
        layer = PhysicalStructure("BadLayer", thickness=-10.)

    # Lower radius above upper
    with pytest.raises(ValueError):
        layer = PhysicalStructure("BadLayer", radius_upper=10., radius_lower=15.)
    
    # Thickness that would produce a bad layer
    with pytest.raises(ValueError):
        layer = PhysicalStructure("BadLayer", radius_upper=10., thickness=15.)
    