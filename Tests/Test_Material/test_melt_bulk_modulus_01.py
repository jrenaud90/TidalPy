"""Melting materials reach the liquid in every modulus together, and the MatPack melts are complete enough to use.

Without a bulk mixing law, the bulk modulus blends linearly into the liquid's across the weakening law's breakdown band,
as the shear modulus does, rather than keeping the solid's until full melt. Every MatPack material that melts over a
range weakens through it, and every MatPack liquid that a melting material melts into can convect (alpha0 > 0).
"""
import pytest

from TidalPy.Material import available_materials, load_material, make_material

PRESSURE = 1.0e9
SURFACE_PRESSURE = 1.0e5


def _state(material, temperature, pressure=PRESSURE):
    return material.calc_state(pressure, temperature, use_melting=True)


def _temperature_at(material, melt_fraction, pressure=PRESSURE):
    solidus, liquidus = material.calc_melting_range(pressure, use_pressure_melting=False)
    return solidus + melt_fraction * (liquidus - solidus)


def test_the_bulk_modulus_blends_across_the_breakdown_band():
    peridotite = load_material("peridotite")
    weakening = peridotite.weakening.parameters
    band_start = weakening["crit_melt_frac"]
    band_end = band_start + weakening["crit_melt_frac_width"]
    liquid = _state(peridotite, _temperature_at(peridotite, 1.0) + 50.0)

    below = _state(peridotite, _temperature_at(peridotite, 0.5 * band_start))
    middle = _state(peridotite, _temperature_at(peridotite, 0.5 * (band_start + band_end)))
    past = _state(peridotite, _temperature_at(peridotite, band_end + 0.05))

    # The solid's below the band, the liquid's past it, and between the two inside it.
    assert below["adiabatic_bulk_modulus"] > 3.0 * liquid["adiabatic_bulk_modulus"]
    assert past["adiabatic_bulk_modulus"] == pytest.approx(liquid["adiabatic_bulk_modulus"], rel=1.0e-2)
    assert past["bulk_modulus"] == pytest.approx(liquid["bulk_modulus"], rel=1.0e-2)
    assert liquid["adiabatic_bulk_modulus"] < middle["adiabatic_bulk_modulus"] < below["adiabatic_bulk_modulus"]


def test_the_bulk_modulus_is_continuous_through_the_band():
    peridotite = load_material("peridotite")
    weakening = peridotite.weakening.parameters
    band_start = weakening["crit_melt_frac"]
    band_end = band_start + weakening["crit_melt_frac_width"]
    step = 1.0e-4
    solid = _state(peridotite, _temperature_at(peridotite, band_start))["adiabatic_bulk_modulus"]
    liquid = _state(peridotite, _temperature_at(peridotite, band_end))["adiabatic_bulk_modulus"]
    # The linear blend's change per step, doubled for the solid and liquid moduli's own drift with temperature; a jump
    # to the liquid would be about (band_end - band_start) / step = 500 times this.
    largest_change = 2.0 * (solid - liquid) * step / (band_end - band_start)
    previous = None
    melt_fraction = band_start - 0.01
    while melt_fraction < band_end + 0.01:
        modulus = _state(peridotite, _temperature_at(peridotite, melt_fraction))["adiabatic_bulk_modulus"]
        if previous is not None:
            assert abs(modulus - previous) < largest_change
        previous = modulus
        melt_fraction += step


def test_without_a_weakening_law_the_bulk_modulus_steps_at_full_melt():
    config = load_material("peridotite").get_config_dict()
    config["melting"].pop("weakening", None)
    plain = make_material(config)
    solid = _state(plain, _temperature_at(plain, 0.0) - 50.0)
    nearly_molten = _state(plain, _temperature_at(plain, 0.99))
    liquid = _state(plain, _temperature_at(plain, 1.0) + 50.0)
    assert nearly_molten["adiabatic_bulk_modulus"] > 3.0 * liquid["adiabatic_bulk_modulus"]
    assert nearly_molten["shear_modulus"] > 0.0
    assert solid["adiabatic_bulk_modulus"] > 3.0 * liquid["adiabatic_bulk_modulus"]


def _range_melting_materials():
    names = []
    for name in available_materials():
        material = load_material(name)
        if material.liquid is None or material.solid is None:
            continue
        solidus, liquidus = material.calc_melting_range(SURFACE_PRESSURE, use_pressure_melting=False)
        if liquidus > solidus:
            names.append(name)
    return names


RANGE_MELTING = _range_melting_materials()


def test_there_are_range_melting_materials():
    assert {"ammonia_water", "iron_sulfide", "peridotite"} <= set(RANGE_MELTING)


@pytest.mark.parametrize("name", RANGE_MELTING)
def test_a_range_melting_material_weakens(name):
    material = load_material(name)
    assert material.weakening is not None
    # Partway through the range, the shear modulus has left the solid's: not rigid until fully molten.
    solidus, liquidus = material.calc_melting_range(SURFACE_PRESSURE, use_pressure_melting=False)
    band_end = material.weakening.parameters["crit_melt_frac"] + material.weakening.parameters["crit_melt_frac_width"]
    temperature = solidus + min(band_end + 0.05, 0.95) * (liquidus - solidus)
    assert _state(material, temperature, SURFACE_PRESSURE)["shear_modulus"] < 1.0


@pytest.mark.parametrize("name", RANGE_MELTING + ["nitrogen_ice"])
def test_a_melting_material_liquid_can_convect(name):
    liquid = load_material(name).liquid
    assert liquid.eos.parameters["thermal_expansion"] > 0.0
