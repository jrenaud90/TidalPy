"""PhysicsLayer viscosity and partial-melt models (attachment, ownership, binary round trip) and EOS model I/O."""
import pytest

from TidalPy.Material.eos import InterpolatedEOS
from TidalPy.Material.eos.material_eos import ConstantDensityEOS
from TidalPy.PartialMelt import make_partial_melt
from TidalPy.Structures.layers.base import BaseLayer
from TidalPy.Structures.layers.gas import GasLayer
from TidalPy.Structures.layers.physics import PhysicsLayer
from TidalPy.Structures.layers.solidliquid import SolidLiquidLayer
from TidalPy.Viscosity import make_viscosity


def _layer_with_material():
    layer = PhysicsLayer("mantle", 0, 0.0, 6.371e6, 4.0e24)
    # The viscosity and partial-melt models belong to the material, so the layer needs one to hand them to.
    layer.set_eos(ConstantDensityEOS(shear_modulus_static=6.0e10, bulk_modulus_static=1.3e11))
    return layer


def _build_layer(layer_class):
    layer = layer_class("mantle", 0, 0.0, 1.0e6, 1.0e22, tidal_scale=0.25)
    layer.set_eos(ConstantDensityEOS(shear_modulus_static=6.0e10, bulk_modulus_static=2.0e11))
    layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity_pas": 3.0e19}))
    layer.set_bulk_viscosity(make_viscosity("reference", {
        "reference_viscosity_pas": 5.0e21, "reference_temperature_k": 1400.0}))
    layer.set_partial_melt(make_partial_melt("henning", {"solidus_k": 1500.0, "liquidus_k": 1900.0}))
    layer.is_static = True
    layer.is_incompressible = True
    return layer


def _roundtrip(layer, layer_class, tmp_path):
    path = str(tmp_path / "layer.tpyb")
    layer.save_binary(path)
    loaded = layer_class("placeholder", 0, 0.0, 1.0, 1.0)
    loaded.load_binary(path)
    return loaded


@pytest.mark.parametrize("attach", [(), ("shear",), ("shear", "bulk", "melt")], ids=["none", "shear_only", "all"])
def test_attach_models_sets_flags(attach):
    """Each set flag turns on only for the model attached; shear and bulk viscosity are independent."""
    layer = _layer_with_material()
    if "shear" in attach:
        layer.set_shear_viscosity(make_viscosity("arrhenius"))
    if "bulk" in attach:
        layer.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity_pas": 1.0e30}))
    if "melt" in attach:
        layer.set_partial_melt(make_partial_melt("henning"))
    assert layer.shear_viscosity_set is ("shear" in attach)
    assert layer.bulk_viscosity_set is ("bulk" in attach)
    assert layer.partial_melt_set is ("melt" in attach)


def test_models_are_move_once():
    """A viscosity or partial-melt model moves into the layer, so it cannot be attached twice."""
    layer = _layer_with_material()
    viscosity = make_viscosity("reference")
    layer.set_shear_viscosity(viscosity)
    with pytest.raises(ValueError):
        layer.set_shear_viscosity(viscosity)
    melt = make_partial_melt("spohn")
    layer.set_partial_melt(melt)
    with pytest.raises(ValueError):
        layer.set_partial_melt(melt)


@pytest.mark.parametrize("layer_class", [PhysicsLayer, SolidLiquidLayer], ids=["physics", "solidliquid"])
def test_strength_models_binary_roundtrip(layer_class, tmp_path):
    """Viscosity and partial-melt models, tidal_scale, and solver flags survive a binary round trip."""
    loaded = _roundtrip(_build_layer(layer_class), layer_class, tmp_path)
    assert loaded.shear_viscosity_set
    assert loaded.bulk_viscosity_set
    assert loaded.partial_melt_set
    assert loaded.tidal_scale == 0.25
    assert loaded.is_static
    assert loaded.is_incompressible
    assert loaded.shear_modulus_static == pytest.approx(6.0e10)


def test_unset_strength_models_stay_unset(tmp_path):
    """A bare layer loads back with no models, no EOS, and an unset tidal scale."""
    loaded = _roundtrip(PhysicsLayer("bare", 0, 0.0, 1.0e6, 1.0e22), PhysicsLayer, tmp_path)
    assert not loaded.shear_viscosity_set
    assert not loaded.bulk_viscosity_set
    assert not loaded.partial_melt_set
    assert not loaded.eos_set
    assert loaded.tidal_scale is None


@pytest.mark.parametrize(
    "layer_class",
    [BaseLayer, PhysicsLayer, SolidLiquidLayer, GasLayer],
    ids=["base", "physics", "solidliquid", "gas"])
def test_eos_model_binary_roundtrip(layer_class, tmp_path):
    """Every layer class saves its material EOS model, including an interpolated model's optional tables."""
    layer = layer_class("mantle", 0, 0.0, 2.0e6, 1.0e22)
    layer.set_eos(InterpolatedEOS(
        radius=[0.0, 1.0e6, 2.0e6],
        density=[5000.0, 4000.0, 3000.0],
        shear_modulus=[1.0e11, 8.0e10, 6.0e10]))
    expected = layer.get_config_dict()["material"]
    loaded = _roundtrip(layer, layer_class, tmp_path)
    assert loaded.eos_set
    assert loaded.get_config_dict()["material"] == expected
