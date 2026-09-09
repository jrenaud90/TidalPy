"""
Tests for the PhysicsLayer viscosity + partial-melt model holders (pipeline
increment A): attaching shear/bulk viscosity and a partial-melt model, the
``*_set`` flags, and move-once ownership semantics.

These models supply the pre-melt viscosities and the melt weakening consumed by
the world EOS solve's per-layer viscoelastic state computation.

Requires the Cython extensions to be compiled first::

    uv pip install -v <repo_root>
"""

import pytest


def _imports():
    from TidalPy.structures_x.layers.physics import PhysicsLayer
    from TidalPy.viscosity_x import make_viscosity
    from TidalPy.partial_melt_x import make_partial_melt
    return PhysicsLayer, make_viscosity, make_partial_melt


def _layer(PhysicsLayer):
    return PhysicsLayer("mantle", 0, 0.0, 6.371e6, 4.0e24,
                        shear_modulus_static_pa=6.0e10, bulk_modulus_static_pa=1.3e11)


def test_flags_start_false():
    PhysicsLayer, _, _ = _imports()
    layer = _layer(PhysicsLayer)
    assert layer.shear_viscosity_set is False
    assert layer.bulk_viscosity_set is False
    assert layer.partial_melt_set is False


def test_attach_models_sets_flags():
    PhysicsLayer, make_viscosity, make_partial_melt = _imports()
    layer = _layer(PhysicsLayer)
    layer.set_shear_viscosity(make_viscosity("arrhenius"))
    layer.set_bulk_viscosity(make_viscosity("constant", {"reference_viscosity": 1.0e30}))
    layer.set_partial_melt(make_partial_melt("henning"))
    assert layer.shear_viscosity_set is True
    assert layer.bulk_viscosity_set is True
    assert layer.partial_melt_set is True


def test_models_are_move_once():
    PhysicsLayer, make_viscosity, make_partial_melt = _imports()
    layer = _layer(PhysicsLayer)
    visc = make_viscosity("reference")
    layer.set_shear_viscosity(visc)
    with pytest.raises(ValueError):
        layer.set_shear_viscosity(visc)   # already moved

    melt = make_partial_melt("spohn")
    layer.set_partial_melt(melt)
    with pytest.raises(ValueError):
        layer.set_partial_melt(melt)


def test_independent_shear_and_bulk_viscosity():
    PhysicsLayer, make_viscosity, _ = _imports()
    layer = _layer(PhysicsLayer)
    layer.set_shear_viscosity(make_viscosity("arrhenius"))
    assert layer.shear_viscosity_set is True
    assert layer.bulk_viscosity_set is False   # bulk independent / still unset
