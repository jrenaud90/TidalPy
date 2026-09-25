"""Material EOS constructors reject parameters with no physical meaning instead of solving to a nonsense planet."""
import math

import pytest

from TidalPy.Material_x.eos import make_material_eos


@pytest.mark.parametrize("model", ("constant", "birch_murnaghan", "vinet"))
@pytest.mark.parametrize("density", (-1000.0, 0.0, math.inf, math.nan))
def test_reference_density_must_be_positive_and_finite(model, density):
    with pytest.raises(ValueError, match="reference density"):
        make_material_eos(model, {"reference_density_kg_m3": density})


@pytest.mark.parametrize("model", ("birch_murnaghan", "vinet"))
@pytest.mark.parametrize("bulk_modulus", (-1.0e11, 0.0, math.inf))
def test_reference_bulk_modulus_must_be_positive_and_finite(model, bulk_modulus):
    with pytest.raises(ValueError, match="reference bulk modulus"):
        make_material_eos(model, {"reference_bulk_modulus_pa": bulk_modulus})


@pytest.mark.parametrize("key", ("shear_modulus_static_pa", "bulk_modulus_static_pa"))
def test_static_moduli_must_not_be_negative(key):
    with pytest.raises(ValueError, match="non-negative"):
        make_material_eos("constant", {key: -1.0})


def test_interpolated_density_must_be_positive():
    with pytest.raises(ValueError, match="density"):
        make_material_eos("interpolate", {"radius_m": [0.0, 1.0e6], "density_kg_m3": [5000.0, -1.0]})
