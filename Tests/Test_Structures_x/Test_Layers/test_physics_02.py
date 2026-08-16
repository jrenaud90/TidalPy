"""Binary round-trip of the PhysicsLayer strength models and classification fields.

Covers the serialized fields beyond the geometry scalars: the attached shear/bulk
viscosity and partial-melt models (each written as a presence flag plus the model's
own recursive record), the ``tidal_scale_method`` selector, and the radial-solver
classification flags (``is_solid``/``is_static``/``is_incompressible``).
"""
import math

import pytest


def _import_layers():
    try:
        from TidalPy.structures_x.layers.physics import PhysicsLayer
        from TidalPy.structures_x.layers.solidliquid import SolidLiquidLayer
        return PhysicsLayer, SolidLiquidLayer
    except ImportError as import_error:
        raise ImportError("structures_x layer extensions are not compiled.") from import_error


def _import_models():
    from TidalPy.viscosity_x.viscosity import make_viscosity
    from TidalPy.partial_melt_x.partial_melt import make_partial_melt
    return make_viscosity, make_partial_melt


def _build_layer(cls):
    make_viscosity, make_partial_melt = _import_models()
    layer = cls(
        "mantle", 0, 0.0, 1.0e6, 1.0e22,
        shear_modulus_static_pa=6.0e10,
        bulk_modulus_static_pa=2.0e11,
        tidal_scale_method="tidal_timescale")
    layer.set_shear_viscosity(make_viscosity("constant", {"reference_viscosity": 3.0e19}))
    layer.set_bulk_viscosity(make_viscosity("reference", {
        "reference_viscosity": 5.0e21, "reference_temperature": 1400.0}))
    layer.set_partial_melt(make_partial_melt("henning", {
        "solidus_k": 1500.0, "liquidus_k": 1900.0}))
    layer.is_static = True
    layer.is_incompressible = True
    return layer


@pytest.mark.parametrize("cls_index", [0, 1], ids=["physics", "solidliquid"])
def test_strength_models_binary_roundtrip(cls_index, tmp_path):
    """Attached viscosity and partial-melt models survive a save/load round trip."""
    classes = _import_layers()
    original = _build_layer(classes[cls_index])
    path = str(tmp_path / "layer.tpyb")
    original.save_binary(path)

    loaded = classes[cls_index]("placeholder", 0, 0.0, 1.0, 1.0)
    loaded.load_binary(path)

    assert loaded.shear_viscosity_set
    assert loaded.bulk_viscosity_set
    assert loaded.partial_melt_set
    assert loaded.tidal_scale_method == original.tidal_scale_method
    assert loaded.is_static
    assert loaded.is_incompressible
    assert loaded.shear_modulus_static == pytest.approx(6.0e10)


def test_unset_strength_models_stay_unset(tmp_path):
    """A layer without strength models loads back with them unset (presence flags off)."""
    PhysicsLayer, _ = _import_layers()
    original = PhysicsLayer("bare", 0, 0.0, 1.0e6, 1.0e22)
    path = str(tmp_path / "bare.tpyb")
    original.save_binary(path)

    loaded = PhysicsLayer("placeholder", 0, 0.0, 1.0, 1.0)
    loaded.load_binary(path)

    assert not loaded.shear_viscosity_set
    assert not loaded.bulk_viscosity_set
    assert not loaded.partial_melt_set
    assert loaded.tidal_scale_method == original.tidal_scale_method
