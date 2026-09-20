"""
Love numbers of a compressible world against the EOS slice count.

The world Love solve reads nothing off the slice grid: the structure comes from the dense EOS solution and the
density and complex moduli from each layer's material state at the exact radius the integrator asks for. The
Love numbers therefore do not depend on slices_per_layer at all, rather than depending on it weakly, and they
converge with the integration tolerance alone.

Requires the Cython extensions to be compiled first::

    uv pip install -v <repo_root>
"""

import cmath
import math

import pytest

from TidalPy.constants import G
from TidalPy.structures_x import build_world

_FREQUENCY = 2.0 * math.pi / 86400.0


def _compressible_world():
    return build_world({
        "schema_version": "0.2.0", "name": "bm", "type": "terrestrial", "radius_m": 6.371e6, "mass_kg": 6.0e24,
        "layers": {
            "core": {"class": "physics", "type": "iron", "layer_index": 0, "radius_fraction": 0.55,
                     "shear_modulus_static_pa": 1.0e11, "bulk_modulus_static_pa": 5.0e11,
                     "eos": {"model": "birch_murnaghan", "reference_density_kg_m3": 8300.0,
                             "reference_bulk_modulus_pa": 1.6e11, "bulk_modulus_derivative": 5.0},
                     "shear_rheology": {"model": "maxwell"},
                     "shear_viscosity": {"model": "constant", "reference_viscosity_pas": 1.0e24}},
            "mantle": {"class": "physics", "type": "mantle_rock", "layer_index": 1, "radius_fraction": 1.0,
                       "shear_modulus_static_pa": 7.0e10, "bulk_modulus_static_pa": 2.0e11,
                       "eos": {"model": "birch_murnaghan", "reference_density_kg_m3": 3300.0,
                               "reference_bulk_modulus_pa": 1.3e11, "bulk_modulus_derivative": 4.0},
                       "shear_rheology": {"model": "maxwell"},
                       "shear_viscosity": {"model": "constant", "reference_viscosity_pas": 1.0e21}}}})


def _k2(world, slices_per_layer, rtol):
    eos = world.solve_eos(G_to_use=G, slices_per_layer=slices_per_layer)
    assert eos["success"] and not eos["max_iters_hit"]
    result = world.solve_love_numbers(frequency=_FREQUENCY, degree_l=2, rtol=rtol, atol=rtol * 1.0e-4)
    assert result["success"], result["message"]
    return world.love_number_k


@pytest.mark.parametrize("slices_per_layer, rel_tol", [(25, 1.0e-7), (50, 1.0e-8), (400, 1.0e-8)])
def test_love_number_of_a_compressible_world_does_not_depend_on_the_slice_count(slices_per_layer, rel_tol):
    """Held from when the density was interpolated between slices; now the agreement is exact."""
    world = _compressible_world()
    reference = _k2(world, 200, 1.0e-11)
    coarse = _k2(world, slices_per_layer, 1.0e-11)
    assert cmath.isclose(coarse, reference, rel_tol=rel_tol), (coarse, reference)


@pytest.mark.parametrize("slices_per_layer", [25, 50, 400])
def test_love_number_of_a_compressible_world_is_bit_identical_across_slice_counts(slices_per_layer):
    """No quantity the radial solve reads is sampled onto the slice grid, so the slice count cannot move k2.

    This is the guarantee the grid-free consumers buy, and it is stronger than the convergence test above:
    not "the interpolation error is small" but "there is no interpolation".
    """
    world = _compressible_world()
    reference = _k2(world, 200, 1.0e-11)
    assert _k2(world, slices_per_layer, 1.0e-11) == reference


def test_love_number_of_a_compressible_world_converges_with_the_tolerance():
    world = _compressible_world()
    reference = _k2(world, 100, 1.0e-11)
    for rtol, expected in ((1.0e-5, 1.0e-6), (1.0e-6, 1.0e-7), (1.0e-8, 1.0e-8)):
        k2 = _k2(world, 100, rtol)
        assert abs(k2.real - reference.real) / abs(reference.real) < expected, (rtol, k2, reference)
        assert abs(k2.imag - reference.imag) / abs(reference) < expected, (rtol, k2, reference)
