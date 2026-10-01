"""BaseWorld and GasGiantWorld: layer ownership, geometry checks, radiogenic heating, binary round-trip."""
import pytest

from TidalPy.Material.eos.material_eos import ConstantDensityEOS
from TidalPy.Radiogenics.radiogenics import FixedRadiogenics
from TidalPy.Structures.layers.base import BaseLayer
from TidalPy.Structures.layers.gas import GasLayer
from TidalPy.Structures.layers.solidliquid import SolidLiquidLayer
from TidalPy.Structures.worlds.base import BaseWorld
from TidalPy.Structures.worlds.gasgiant import GasGiantWorld


# Earth-like two-layer geometry (MKS).
_R_CMB   = 3.485e6
_R_SURF  = 6.371e6
_M_CORE  = 1.932e24
_M_MANT  = 4.040e24
_M_TOT   = _M_CORE + _M_MANT

# Positional layer arguments: name, index, radius_inner, radius_outer, mass.
_CORE_ARGS   = ("core", 0, 0.0, _R_CMB, _M_CORE)
_MANTLE_ARGS = ("mantle", 1, _R_CMB, _R_SURF, _M_MANT)
_ENVELOPE_ARGS = ("envelope", 0, 0.0, 7.0e7, 1.898e27)


def _make_two_layer_world():
    world = BaseWorld("Earth", _R_SURF, _M_TOT, world_type="terrestrial")
    world.add_layer(SolidLiquidLayer(*_CORE_ARGS, material_name="iron"))
    world.add_layer(SolidLiquidLayer(*_MANTLE_ARGS, material_name="perovskite"))
    return world


def test_add_layers_and_count():
    """Two added layers are counted and validate."""
    world = _make_two_layer_world()
    assert world.num_layers == 2
    assert world.validate_layers() is True


def test_add_layer_consumes_wrapper():
    """A layer wrapper cannot be added twice."""
    world = BaseWorld("W", _R_SURF, _M_TOT)
    layer = SolidLiquidLayer(*_CORE_ARGS)
    world.add_layer(layer)
    with pytest.raises(ValueError):
        world.add_layer(layer)


def test_added_layer_wrapper_stays_usable():
    """The wrapper handed to add_layer becomes a view of the world-owned layer."""
    world = BaseWorld("W", _R_CMB, _M_CORE)
    layer = SolidLiquidLayer(*_CORE_ARGS)
    layer.set_eos(ConstantDensityEOS(reference_density=9000.0))
    world.add_layer(layer)

    assert layer.name == "core"
    assert layer.radius_outer == _R_CMB
    layer.is_static = False
    assert world.core.is_static is False
    assert world.solve_eos()["success"]
    assert layer.get_density(0.5 * _R_CMB) == pytest.approx(9000.0)
    assert layer.get_density(0.5 * _R_CMB) == world.core.get_density(0.5 * _R_CMB)


def test_attached_model_wrappers_raise_instead_of_reading_a_moved_object():
    """Model wrappers whose C++ object moved into a layer raise on use."""
    from TidalPy.Rheology import Maxwell
    from TidalPy.Viscosity import make_viscosity
    from TidalPy.PartialMelt import make_partial_melt

    layer = SolidLiquidLayer("mantle", 0, 0.0, _R_SURF, _M_TOT)
    layer.set_eos(ConstantDensityEOS(reference_density=5000.0))

    # Spec-driven models (rheology, viscosity) are copied in, so their wrappers stay usable.
    rheology = Maxwell()
    layer.set_shear_rheology(rheology)
    assert rheology.model_name == "maxwell"

    # A spec-driven model (viscosity) is copied in, so its wrapper stays usable.
    viscosity = make_viscosity("constant", {"reference_viscosity_pas": 1.0e20})
    layer.set_shear_viscosity(viscosity)
    assert viscosity.calc_viscosity(1500.0, 1.0e9) == pytest.approx(1.0e20)

    partial_melt = make_partial_melt("henning")
    layer.set_partial_melt(partial_melt)
    with pytest.raises(RuntimeError, match="took ownership"):
        partial_melt.calc_melt_fraction(1700.0)


def test_add_layer_discontinuity_raises():
    """A layer leaving a radial gap is rejected but not consumed."""
    world = BaseWorld("W", _R_SURF, _M_TOT)
    # The first layer must start at radius 0.
    bad = SolidLiquidLayer("mantle", 0, _R_CMB, _R_SURF, _M_MANT)
    with pytest.raises(ValueError):
        world.add_layer(bad)
    world2 = BaseWorld("W2", _R_SURF, _M_MANT)
    world2.add_layer(SolidLiquidLayer(*_CORE_ARGS))
    world2.add_layer(bad)
    assert world2.num_layers == 2


def test_calc_total_mass():
    """calc_total_mass sums the layer masses."""
    world = _make_two_layer_world()
    assert world.calc_total_mass() == pytest.approx(_M_TOT, rel=1e-9)


def test_mixed_layer_types():
    """A world accepts layers of different subclasses."""
    world = BaseWorld("Mixed", _R_SURF, _M_TOT)
    world.add_layer(BaseLayer(*_CORE_ARGS))
    world.add_layer(BaseLayer(*_MANTLE_ARGS))
    assert world.num_layers == 2
    assert world.calc_total_mass() == pytest.approx(_M_TOT, rel=1e-9)


def test_internal_heating_zero_without_radiogenics():
    """Internal heating is zero when no layer has radiogenics."""
    world = _make_two_layer_world()
    assert world.calc_internal_heating(0.0) == pytest.approx(0.0)


def test_internal_heating_with_radiogenics():
    """Internal heating is the mantle's radiogenic rate times its mass."""
    world = BaseWorld("Earth", _R_SURF, _M_TOT)
    mantle = SolidLiquidLayer(*_MANTLE_ARGS)
    mantle.set_radiogenics(FixedRadiogenics(fixed_heat_production=1.0e-11))
    world.add_layer(SolidLiquidLayer(*_CORE_ARGS))
    world.add_layer(mantle)
    assert world.calc_internal_heating(0.0) == pytest.approx(1.0e-11 * _M_MANT, rel=1e-9)


def test_layered_world_binary_roundtrip(tmp_path):
    """Binary save/load restores the world, its layers, and their sub-models."""
    world = BaseWorld(
        "Earth",
        _R_SURF,
        _M_TOT,
        world_type="terrestrial",
        albedo=0.31,
        obliquity=0.41,
    )
    mantle = SolidLiquidLayer(*_MANTLE_ARGS, material_name="perovskite")
    mantle.set_radiogenics(FixedRadiogenics(fixed_heat_production=2.0e-11))
    world.add_layer(SolidLiquidLayer(*_CORE_ARGS, material_name="iron"))
    world.add_layer(mantle)
    heating_before = world.calc_internal_heating(0.0)

    path = str(tmp_path / "world.tpyb")
    world.save_binary(path)
    loaded = BaseWorld("placeholder", 1.0, 1.0)
    loaded.load_binary(path)

    assert loaded.name        == "Earth"
    assert loaded.world_type  == "terrestrial"
    assert loaded.albedo      == pytest.approx(0.31)
    assert loaded.obliquity   == pytest.approx(0.41)
    assert loaded.num_layers  == 2
    assert loaded.calc_total_mass() == pytest.approx(_M_TOT, rel=1e-9)
    assert loaded.calc_internal_heating(0.0) == pytest.approx(heating_before, rel=1e-12)
    cfg = loaded.get_config_dict()
    assert list(cfg["layers"]) == ["core", "mantle"]
    assert cfg["layers"]["mantle"]["radius_outer_m"] == pytest.approx(_R_SURF)
    assert "radius_inner_m" not in cfg["layers"]["mantle"]
    assert "radiogenics" in cfg["layers"]["mantle"]


def test_gasgiant_construction_and_type():
    """A GasGiantWorld reports its type and accepts a gas layer."""
    gas_giant = GasGiantWorld("Jupiter", 7.0e7, 1.898e27)
    assert gas_giant.world_type == "gasgiant"
    gas_giant.add_layer(GasLayer(*_ENVELOPE_ARGS))
    assert gas_giant.num_layers == 1


def test_gasgiant_binary_roundtrip(tmp_path):
    """Binary save/load restores a GasGiantWorld."""
    gas_giant = GasGiantWorld("Jupiter", 7.0e7, 1.898e27)
    gas_giant.add_layer(GasLayer(*_ENVELOPE_ARGS))
    path = str(tmp_path / "gasgiant.tpyb")
    gas_giant.save_binary(path)
    loaded = GasGiantWorld("placeholder", 1.0, 1.0)
    loaded.load_binary(path)
    assert loaded.name       == "Jupiter"
    assert loaded.world_type == "gasgiant"
    assert loaded.num_layers == 1


@pytest.mark.parametrize(
    "make_world, parent_class",
    [
        (_make_two_layer_world, BaseWorld),
        (lambda: GasGiantWorld("Jupiter", 7.0e7, 1.898e27), BaseWorld),
    ],
    ids=["layered_is_base", "gasgiant_is_layered"],
)
def test_world_class_hierarchy(make_world, parent_class):
    """World classes subclass their parent world class."""
    assert isinstance(make_world(), parent_class)
