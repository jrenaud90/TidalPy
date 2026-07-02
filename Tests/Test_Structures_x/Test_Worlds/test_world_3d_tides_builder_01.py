"""Builder wiring of the 2D tidal potential model onto a world's 3D dissipation path.

The world builder attaches a tidal potential model (for the on-demand 3D stress/strain/heating
path) when a layered world runs the rheology dissipation model. These tests check that the
``[tides] tidal_potential_model`` key selects the model, that the default is used when omitted,
and that the potential is not attached for the analytic dissipation models (which have no
depth-resolved solution) nor for a layerless star.
"""
import numpy as np

from TidalPy.constants import G
from TidalPy.structures_x import build_world
from TidalPy.structures_x.worlds.layered import LayeredWorld


# Io-Jupiter-like orbital state.
_HOST_MASS = 1.898e27
_SMA = 4.2e8
_N = 2.05e-5
_SPIN = 1.5 * _N
_ECC = 0.05


def _rheology_config(tides_extra=None):
    tides = {
        "global_tidal_model": "rheology", "max_degree_l": 2,
        "eccentricity_trunc_lvl": 2, "obliquity_trunc_lvl": "off",
    }
    if tides_extra:
        tides.update(tides_extra)
    return {
        "name": "test_terr", "type": "terrestrial", "radius_m": 1.0e6, "mass_kg": 2.0e22,
        "tides": tides,
        "layers": {
            "mantle": {
                "class": "physics", "type": "mantle_rock", "radius_fraction": 1.0,
                "shear_modulus_static_pa": 5.0e10, "bulk_modulus_static_pa": 1.0e11,
                "shear_viscosity": {"model": "constant", "reference_viscosity": 1.0e19},
                "bulk_viscosity": {"model": "constant", "reference_viscosity": 1.0e19},
                "shear_rheology": {"model": "maxwell"},
                "bulk_rheology": {"model": "elastic"},
            }
        },
    }


def test_default_potential_model_attached_for_rheology():
    world = build_world(_rheology_config())
    assert isinstance(world, LayeredWorld)
    assert world.tidal_potential_model_set


def test_selected_potential_model_runs_3d_path():
    world = build_world(_rheology_config({"tidal_potential_model": "nsr_modes"}))
    assert world.tidal_potential_model_set
    world.solve_eos(G_to_use=G)
    heat = world.get_3d_tidal_heating(_N, _SPIN, _ECC, 0.0, _SMA, _HOST_MASS, 0.5e6, 1.0, 0.5, 0.0)
    assert np.isfinite(heat)
    assert heat > 0.0


def test_obliquity_potential_model_selectable():
    world = build_world(_rheology_config({"tidal_potential_model": "nsr_modes_med_obliquity"}))
    assert world.tidal_potential_model_set
    world.solve_eos(G_to_use=G)
    # Nonzero obliquity activates the m = 1 modes; the 3D path must still evaluate finitely.
    heat = world.get_3d_tidal_heating(_N, _SPIN, _ECC, 0.2, _SMA, _HOST_MASS, 0.5e6, 1.0, 0.5, 0.0)
    assert np.isfinite(heat)


def test_analytic_model_gets_no_potential():
    """cpl/ctl/ctl_q have no depth-resolved solution, so no potential model is attached."""
    cfg = _rheology_config()
    cfg["tides"]["global_tidal_model"] = "cpl"
    cfg["tides"]["fixed_k"] = [0.3]
    cfg["tides"]["fixed_q"] = [50.0]
    world = build_world(cfg)
    assert not world.tidal_potential_model_set


def test_star_gets_no_potential_model():
    """A layerless star exposes no set_tidal_potential_model; the builder must skip it cleanly."""
    star_cfg = {
        "name": "test_star", "type": "star", "radius_m": 7.0e8, "mass_kg": 2.0e30,
        "tides": {"global_tidal_model": "fixed_q", "fixed_k": [0.03], "fixed_q": [1.0e6]},
    }
    star = build_world(star_cfg)
    assert not hasattr(star, "tidal_potential_model_set") or not star.tidal_potential_model_set
