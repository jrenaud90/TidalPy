"""The point-wise strain, stress, and heating kernel against the frozen classic results."""
from pathlib import Path

import numpy as np

from TidalPy.Tides.multilayer.stress_strain import strain_stress_heating_point

# Classic degree-2 `calculate_strain_stress` and heating bilinear form at 300 inputs drawn with default_rng(99),
# from the 0.8.0 snapshot 8b8e0b12 (solid compressible, the case the classic code supports).
FROZEN_PATH = Path(__file__).parent / "frozen" / "test_point_kernel_01.npz"
with np.load(FROZEN_PATH, allow_pickle=False) as frozen_file:
    LEGACY_REFERENCE = {key: frozen_file[key] for key in frozen_file.files}


def test_strain_stress_heating_matches_legacy():
    """Strains, stresses, and heating match the classic results to machine precision."""
    worst_e = worst_s = worst_h = 0.0
    for sample_i in range(LEGACY_REFERENCE["radii"].size):
        y = LEGACY_REFERENCE["radial_functions"][sample_i]
        shear = complex(LEGACY_REFERENCE["shear_moduli"][sample_i])
        bulk = complex(LEGACY_REFERENCE["bulk_moduli"][sample_i])
        radius = float(LEGACY_REFERENCE["radii"][sample_i])
        pot6 = tuple(LEGACY_REFERENCE["potential_terms"][sample_i])
        theta = float(LEGACY_REFERENCE["colatitudes"][sample_i])
        frequency = float(LEGACY_REFERENCE["frequencies"][sample_i])

        e_ref = LEGACY_REFERENCE["classic_strains"][sample_i]
        s_ref = LEGACY_REFERENCE["classic_stresses"][sample_i]
        # The classic value is the bilinear form alone; the kernel returns |frequency| / 2 times it.
        h_ref = float(LEGACY_REFERENCE["classic_heating_bilinear"][sample_i])
        h_ref *= 0.5 * abs(frequency)
        e_c, s_c, h_c = strain_stress_heating_point(
            np.ascontiguousarray(y, dtype=np.complex128),
            shear,
            bulk,
            radius,
            2.0,
            frequency,
            True,
            False,
            pot6,
            theta)
        worst_e = max(worst_e, float(np.max(np.abs(e_c - e_ref) / (np.abs(e_ref) + 1e-30))))
        worst_s = max(worst_s, float(np.max(np.abs(s_c - s_ref) / (np.abs(s_ref) + 1e-30))))
        worst_h = max(worst_h, abs(h_c - h_ref) / (abs(h_ref) + 1e-30))
    assert worst_e < 1e-11, f"strain differs: {worst_e:.3e}"
    assert worst_s < 1e-11, f"stress differs: {worst_s:.3e}"
    assert worst_h < 1e-10, f"heating differs: {worst_h:.3e}"


def test_liquid_returns_nan():
    """The shear kernel is solid-only, so a liquid point returns NaN."""
    y = np.ones(6, dtype=np.complex128)
    e, s, h = strain_stress_heating_point(
        y,
        0.0 + 0j,
        1.3e11 + 0j,
        1e6,
        2.0,
        1.0e-5,
        False,
        False,
        (1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
        1.0)
    assert np.all(np.isnan(e)) and np.isnan(h)
