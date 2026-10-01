"""A molten band thinner than the EOS reporting grid. The first melt of a heated mantle forms where its conducting
top boundary layer meets the adiabat, and it can lie wholly between two slices; it is found, and the Love solve
treats it as a static liquid, whatever slices_per_layer is."""
import pytest

from TidalPy.Structures import build_world

# Monteux et al. (2016) peridotite melting curves, and mantle-silicate Anderson-Gruneisen constants.
_MONTEUX = {
    "solidus_k": 1661.2, "solidus_simon_a_pa": 1.336e9, "solidus_simon_c": 7.437,
    "solidus_transition_pressure_pa": 20.0e9, "solidus_high_k": 2081.8, "solidus_high_simon_a_pa": 1.0169e11,
    "solidus_high_simon_c": 1.226,
    "liquidus_k": 1982.1, "liquidus_simon_a_pa": 6.594e9, "liquidus_simon_c": 5.374,
    "liquidus_transition_pressure_pa": 20.0e9, "liquidus_high_k": 2006.8, "liquidus_high_simon_a_pa": 3.465e10,
    "liquidus_high_simon_c": 1.844,
}
_SURFACE_TEMPERATURE = 815.0   # [K] an Earth's equilibrium temperature at 0.1 AU from the Sun
_FREQUENCY = 6.3e-6            # [rad s-1] about the orbital motion there
_COARSE_SLICES = 100
_FINE_SLICES = 3000


def _earth(mantle_temperature, slices_per_layer):
    config = build_world("earth_simple").get_config_dict()
    config.setdefault("eos_solver", {})["slices_per_layer"] = slices_per_layer
    mantle = config["layers"]["mantle"]
    mantle["temperature_k"] = mantle_temperature
    mantle["material"]["anderson_gruneisen_parameter"] = 5.5
    mantle["material"]["anderson_gruneisen_exponent"] = 1.4
    mantle["material"]["partial_melt"].update(_MONTEUX)
    world = build_world(config)
    result = world.solve_eos(solve_temperature=True, surface_temperature=_SURFACE_TEMPERATURE)
    assert result["success"], result["message"]
    return world


def _upper_mantle_bands(world):
    """The molten stretches in the upper half of the mantle, as (radius_inner, radius_outer) [m]."""
    middle = 0.5 * (world.mantle.radius_inner + world.mantle.radius_outer)
    return [(inner, outer) for name, inner, outer in world.molten_regions if (name == "mantle") and (inner > middle)]


@pytest.mark.parametrize("mantle_temperature", [1875.0, 1890.0])
def test_a_band_thinner_than_the_slice_spacing_is_found(mantle_temperature):
    """The near-surface band is a few km thick, well under the 29 km spacing of 100 slices in this mantle, and the
    coarse grid finds the same band as a grid 30 times finer."""
    coarse = _earth(mantle_temperature, _COARSE_SLICES)
    fine = _earth(mantle_temperature, _FINE_SLICES)
    coarse_bands = _upper_mantle_bands(coarse)
    fine_bands = _upper_mantle_bands(fine)
    assert len(fine_bands) == 1
    assert len(coarse_bands) == 1
    slice_spacing = (coarse.mantle.radius_outer - coarse.mantle.radius_inner) / (_COARSE_SLICES - 1)
    band_inner, band_outer = coarse_bands[0]
    assert 0.0 < band_outer - band_inner < 0.5 * slice_spacing
    assert band_inner == pytest.approx(fine_bands[0][0], abs=1.0)
    assert band_outer == pytest.approx(fine_bands[0][1], abs=1.0)


@pytest.mark.parametrize("mantle_temperature", [1875.0, 1890.0])
def test_the_love_solve_does_not_depend_on_the_slice_count(mantle_temperature):
    """With the band missed, the shooting method either integrated it as a solid at the liquid-shear floor and failed,
    or returned a Love number for a mantle without the band (22 times the tidal heating at 1875 K). Both grids now
    give the same k2, to the radial solve's own tolerance (rtol 1e-6 by default); the band edges agree to a mm."""
    love_numbers = []
    for slices in (_COARSE_SLICES, _FINE_SLICES):
        world = _earth(mantle_temperature, slices)
        world.solve_love_numbers(frequency=_FREQUENCY, degree_l=2)
        assert world.love_success, world.love_message
        love_numbers.append(complex(world.love_number_k))
    assert love_numbers[0] == pytest.approx(love_numbers[1], rel=1.0e-5)
