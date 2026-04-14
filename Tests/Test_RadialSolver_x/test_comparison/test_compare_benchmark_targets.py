import json
from pathlib import Path

import numpy as np
import pytest

from TidalPy.RadialSolver_x import build_rs_input_homogeneous_layers
from TidalPy.RadialSolver_x.solver import radial_solver
from TidalPy.rheology.models import Andrade, Elastic, Maxwell

TARGETS_PATH = Path(__file__).parent / "data" / "general_benchmarks_targets.json"

BASE_KWARGS = dict(
    surface_pressure=0,
    degree_l=2,
    solve_for=None,
    core_model=0,
    use_kamata=False,
    starting_radius=0.0,
    start_radius_tolerance=1.0e-5,
    integration_method="DOP853",
    integration_rtol=1.0e-4,
    integration_atol=1.0e-6,
    scale_rtols_bylayer_type=False,
    max_num_steps=500_000,
    expected_size=200,
    max_ram_MB=500,
    max_step=0,
    nondimensionalize=True,
    use_prop_matrix=False,
    verbose=False,
    warnings=False,
    raise_on_fail=False,
    eos_method_bylayer=None,
    eos_integration_method="RK45",
    eos_rtol=1.0e-4,
    eos_atol=1.0e-6,
    eos_pressure_tol=1.0e-2,
    eos_max_iters=50,
    perform_checks=False,
)


def _load_targets():
    with TARGETS_PATH.open("r", encoding="utf-8") as f:
        payload = json.load(f)
    return {item["case"]: item for item in payload["targets"]}


def _complex_array(pairs):
    arr = np.asarray(pairs, dtype=np.float64)
    return arr[..., 0] + 1.0j * arr[..., 1]


def _build_case(case_name):
    kwargs = dict(BASE_KWARGS)

    if case_name == "1layer":
        rs_input = build_rs_input_homogeneous_layers(
            planet_radius=6000.0e3,
            forcing_frequency=np.pi * 2.0 / (86400.0 * 7.5),
            density_tuple=(5400.0,),
            static_bulk_modulus_tuple=(1.0e11,),
            static_shear_modulus_tuple=(50.0e9,),
            bulk_viscosity_tuple=(1.0e18,),
            shear_viscosity_tuple=(1.0e18,),
            layer_type_tuple=("solid",),
            layer_is_static_tuple=(False,),
            layer_is_incompressible_tuple=(False,),
            shear_rheology_model_tuple=(Andrade(),),
            bulk_rheology_model_tuple=(Elastic(),),
            thickness_fraction_tuple=(1.0,),
            volume_fraction_tuple=None,
            slices_tuple=None,
            slice_per_layer=10,
            perform_checks=False,
        )
    elif case_name == "2layer":
        rs_input = build_rs_input_homogeneous_layers(
            planet_radius=6000.0e3,
            forcing_frequency=np.pi * 2.0 / (86400.0 * 0.3),
            density_tuple=(11000.0, 3400.0),
            static_bulk_modulus_tuple=(5.0e11, 1.0e11),
            static_shear_modulus_tuple=(0.0, 50.0e9),
            bulk_viscosity_tuple=(1000.0, 1.0e18),
            shear_viscosity_tuple=(1000.0, 1.0e18),
            layer_type_tuple=("liquid", "solid"),
            layer_is_static_tuple=(False, False),
            layer_is_incompressible_tuple=(False, False),
            shear_rheology_model_tuple=(Elastic(), Andrade()),
            bulk_rheology_model_tuple=(Elastic(), Elastic()),
            thickness_fraction_tuple=(0.4, 0.6),
            volume_fraction_tuple=None,
            slices_tuple=None,
            slice_per_layer=30,
            perform_checks=False,
        )
    elif case_name == "3layer":
        rs_input = build_rs_input_homogeneous_layers(
            planet_radius=6000.0e3,
            forcing_frequency=np.pi * 2.0 / (86400.0 * 0.1),
            density_tuple=(9600.0, 8000.0, 3400.0),
            static_bulk_modulus_tuple=(10.0e11, 5.0e11, 1.0e11),
            static_shear_modulus_tuple=(150.0e9, 0.0, 50.0e9),
            bulk_viscosity_tuple=(1.0e27, 1000.0, 1.0e18),
            shear_viscosity_tuple=(1.0e27, 1000.0, 1.0e18),
            layer_type_tuple=("solid", "liquid", "solid"),
            layer_is_static_tuple=(False, False, False),
            layer_is_incompressible_tuple=(False, False, False),
            shear_rheology_model_tuple=(Andrade(), Elastic(), Andrade()),
            bulk_rheology_model_tuple=(Elastic(), Elastic(), Elastic()),
            thickness_fraction_tuple=(0.15, 0.30, 0.55),
            volume_fraction_tuple=None,
            slices_tuple=None,
            slice_per_layer=30,
            perform_checks=False,
        )
    elif case_name == "4layer":
        rs_input = build_rs_input_homogeneous_layers(
            planet_radius=2574765.0,
            forcing_frequency=np.pi * 2.0 / (86400.0 * 5.0),
            density_tuple=(2565.0, 1250.0, 1122.0, 950.0),
            static_bulk_modulus_tuple=(100.0e9, 25.0e9, 3.10e9, 9.70e9),
            static_shear_modulus_tuple=(50.0e9, 4.0 * 3.24e9, 0.0, 3.24e9),
            bulk_viscosity_tuple=(1.0e27, 2.0 * 3.24e9, 1000.0, 1.00e12),
            shear_viscosity_tuple=(1.0e27, 2.0 * 3.24e9, 1000.0, 1.00e12),
            layer_type_tuple=("solid", "solid", "liquid", "solid"),
            layer_is_static_tuple=(False, False, False, False),
            layer_is_incompressible_tuple=(False, False, False, False),
            shear_rheology_model_tuple=(Maxwell(), Maxwell(), Elastic(), Maxwell()),
            bulk_rheology_model_tuple=(Elastic(), Elastic(), Elastic(), Elastic()),
            thickness_fraction_tuple=(0.80901, 0.05049, 0.04350, 0.09700),
            volume_fraction_tuple=None,
            slices_tuple=None,
            slice_per_layer=20,
            perform_checks=False,
        )

        kwargs["integration_rtol"] = 1.0e-18
        kwargs["integration_atol"] = 1.0e-21
        kwargs["starting_radius"] = 0.0
    else:
        raise ValueError(f"Unknown benchmark case: {case_name}")

    return rs_input, kwargs


@pytest.mark.parametrize("case_name", ("1layer", "2layer", "3layer", "4layer"))
def test_benchmark_targets(case_name):
    targets = _load_targets()[case_name]
    rs_input, kwargs = _build_case(case_name)

    solution = radial_solver(*rs_input, **kwargs)

    assert solution.success, solution.message

    expected_steps = np.asarray(targets["steps_required"], dtype=np.int64)
    expected_love = _complex_array(targets["love"])

    np.testing.assert_array_equal(solution.steps_taken, expected_steps)
    np.testing.assert_allclose(solution.love, expected_love, rtol=1.0e-7, atol=1.0e-10)
