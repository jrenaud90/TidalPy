import TidalPy


from TidalPy.structures import build_world

# Build basic layered io
io_base = build_world('io_simple')


def test_basic_set_temperature():
    # Set the temperature of the Mantle to something where there would be no partial melt
    io_base.Mantle.set_state(temperature=1500.)

    assert io_base.Mantle.temperature == 1500.
    assert io_base.Mantle.rheology.partial_melting_model.melt_fraction == 0.
    assert io_base.Mantle.rheology.melt_fraction is io_base.Mantle.rheology.partial_melting_model.melt_fraction
    assert io_base.Mantle.melt_fraction is io_base.Mantle.rheology.partial_melting_model.melt_fraction

    # Make sure the layer below the mantle is picking up on the changes
    assert io_base.Core.surface_temperature is io_base.Mantle.temperature

    # Make sure the strength information is passing around where it should be
    assert io_base.Mantle.rheology.postmelt_viscosity is \
           io_base.Mantle.rheology.partial_melting_model.postmelt_viscosity
    assert io_base.Mantle.viscosity is io_base.Mantle.rheology.partial_melting_model.postmelt_viscosity
    assert io_base.Mantle.rheology.postmelt_shear_modulus is \
           io_base.Mantle.rheology.partial_melting_model.postmelt_shear_modulus
    assert io_base.Mantle.shear_modulus is io_base.Mantle.rheology.partial_melting_model.postmelt_shear_modulus
    assert io_base.Mantle.rheology.postmelt_compliance is \
           io_base.Mantle.rheology.partial_melting_model.postmelt_compliance
    assert io_base.Mantle.compliance is io_base.Mantle.rheology.partial_melting_model.postmelt_compliance

    premelt_shear = io_base.Mantle.shear_modulus
    premelt_visco = io_base.Mantle.viscosity

    # Try with some melt now
    io_base.Mantle.set_state(temperature=1800.)

    assert io_base.Mantle.melt_fraction > 0.
    assert io_base.Mantle.melt_fraction < 1.
    assert io_base.Mantle.shear_modulus < premelt_shear
    assert io_base.Mantle.viscosity < premelt_visco
