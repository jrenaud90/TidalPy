"""A molten band thinner than the EOS reporting grid. A heated mantle on pressure-dependent melting curves is most
molten at the top of its convecting interior, where the layer's temperature applies at the lowest pressure. Just past
the rheological transition there the interior is a magma ocean (the convection model takes its viscosity at that
point), its conducting lid is far thinner than a slice, and the molten band from that depth up to the lid is only a
few km thick; it is found, and the Love solve treats it as a static liquid, whatever slices_per_layer is."""
import pytest

from TidalPy.Structures import build_world

_SURFACE_TEMPERATURE = 815.0   # [K] an Earth's equilibrium temperature at 0.1 AU from the Sun
_FREQUENCY = 6.3e-6            # [rad s-1] about the orbital motion there
_COARSE_SLICES = 100
_FINE_SLICES = 3000


# Monteux et al. (2016) peridotite melting curves, and mantle-silicate Anderson-Gruneisen constants.
_MONTEUX_SOLIDUS = {
    "model": "simon_glatzel_2", "temperature_k": 1661.2, "simon_a_pa": 1.336e9, "simon_c": 7.437,
    "transition_pressure_pa": 20.0e9, "high_temperature_k": 2081.8, "high_simon_a_pa": 1.0169e11,
    "high_simon_c": 1.226}
_MONTEUX_LIQUIDUS = {
    "model": "simon_glatzel_2", "temperature_k": 1982.1, "simon_a_pa": 6.594e9, "simon_c": 5.374,
    "transition_pressure_pa": 20.0e9, "high_temperature_k": 2006.8, "high_simon_a_pa": 3.465e10,
    "high_simon_c": 1.844}
_DELTA = 5.5
_KAPPA = 1.4
# Silicate and iron thermal constants: conductivity [W m-1 K-1], heat capacity [J kg-1 K-1], expansivity [1/K].
_ROCK_THERMAL = {"thermal_conductivity_w_mk": 3.75, "heat_capacity_j_kgk": 1200.0}
_ROCK_EXPANSION = 5.2e-5
_IRON_THERMAL = {"thermal_conductivity_w_mk": 7.95, "heat_capacity_j_kgk": 840.0}
_IRON_EXPANSION = 1.2e-5


def _thermal_earth_config(pressure_dependent):
    """Bundled earth_simple, ready for a thermal solve: iron cores, and a convecting mantle with silicate thermal
    constants that melts into a Murnaghan melt of 0.2 Pa s through Henning weakening. Its melting curves are a constant
    1600 K solidus and 2000 K liquidus, or with pressure_dependent the Monteux curves, which follow the pressure, and an
    Anderson-Gruneisen expansivity."""
    config = build_world("earth_simple").get_config_dict()
    for name in ("inner_core", "outer_core"):
        layer = config["layers"][name]
        for phase in layer["material"].values():
            if isinstance(phase, dict):
                phase.update(_IRON_THERMAL)
                phase["eos"]["thermal_expansion_1_k"] = _IRON_EXPANSION
        layer["cooling"] = {"model": "off"}
        layer["radiogenics"] = {"model": "off"}
    mantle = config["layers"]["mantle"]
    material = mantle["material"]
    material["solid"].update(_ROCK_THERMAL)
    material["solid"]["eos"]["thermal_expansion_1_k"] = _ROCK_EXPANSION
    material["liquid"] = {
        **_ROCK_THERMAL,
        "eos": {"model": "murnaghan", "reference_density_kg_m3": 2750.0, "reference_bulk_modulus_pa": 2.0e10,
                "bulk_modulus_derivative": 5.0, "thermal_expansion_1_k": _ROCK_EXPANSION},
        "shear_viscosity": {"model": "constant", "reference_viscosity_pas": 0.2}}
    material["melting"] = {"solidus": {"model": "constant", "temperature_k": 1600.0},
                           "liquidus": {"model": "constant", "temperature_k": 2000.0},
                           "weakening": {"model": "henning"}}
    mantle["use_melting"] = True
    mantle["cooling"] = {"model": "convection", "convection_alpha": 1.0, "convection_beta": 1.0 / 3.0,
                         "critical_rayleigh": 1100.0}
    mantle["radiogenics"] = {"model": "isotope", "isotopes": "modern_day_chondritic"}
    if pressure_dependent:
        material["solid"]["eos"]["anderson_gruneisen_parameter"] = _DELTA
        material["solid"]["eos"]["anderson_gruneisen_exponent"] = _KAPPA
        material["melting"]["solidus"] = dict(_MONTEUX_SOLIDUS)
        material["melting"]["liquidus"] = dict(_MONTEUX_LIQUIDUS)
        mantle["use_pressure_melting"] = True
    return config


def _earth(mantle_temperature, slices_per_layer):
    config = _thermal_earth_config(pressure_dependent=True)
    config.setdefault("eos_solver", {})["slices_per_layer"] = slices_per_layer
    config["layers"]["mantle"]["temperature_k"] = mantle_temperature
    world = build_world(config)
    result = world.solve_eos(solve_temperature=True, surface_temperature=_SURFACE_TEMPERATURE)
    assert result["success"], result["message"]
    return world


def _upper_mantle_bands(world):
    """The molten stretches in the upper half of the mantle, as (radius_inner, radius_outer) [m]."""
    middle = 0.5 * (world.mantle.radius_inner + world.mantle.radius_outer)
    return [(inner, outer) for name, inner, outer in world.molten_regions if (name == "mantle") and (inner > middle)]


# Just past the rheological transition (Henning's critical melt fraction, 0.5) at the top of the interior: 0.59 melt at
# 1850 K. Below it (1825 K, 0.49) the interior is a partially molten solid and no stretch is molten.
_BAND_TEMPERATURES = (1850.0, 1855.0)


@pytest.mark.parametrize("mantle_temperature", _BAND_TEMPERATURES)
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


@pytest.mark.parametrize("mantle_temperature", _BAND_TEMPERATURES)
def test_the_love_solve_does_not_depend_on_the_slice_count(mantle_temperature):
    """A missed band would either be integrated as a solid of near-zero rigidity, failing the shooting method, or
    leave a Love number for a mantle without the band. Both grids give the same k2, to the radial solve's own
    tolerance (rtol 1e-6 by default)."""
    love_numbers = []
    for slices in (_COARSE_SLICES, _FINE_SLICES):
        world = _earth(mantle_temperature, slices)
        world.solve_love_numbers(frequency=_FREQUENCY, degree_l=2)
        assert world.love_success, world.love_message
        love_numbers.append(complex(world.love_number_k))
    assert love_numbers[0] == pytest.approx(love_numbers[1], rel=1.0e-5)
