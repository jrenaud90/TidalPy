"""The point-wise strain, stress, and heating kernel against the frozen classic results."""
from pathlib import Path

import numpy as np

from TidalPy.Tides.multilayer.stress_strain import strain_stress_heating_point, volumetric_heating

# Classic degree-2 `calculate_strain_stress` and heating bilinear form at 300 inputs drawn with default_rng(99),
# from the 0.8.0 snapshot 8b8e0b12 (solid compressible, the case the classic code supports).
FROZEN_PATH = Path(__file__).parent / "frozen" / "test_point_kernel_01.npz"
with np.load(FROZEN_PATH, allow_pickle=False) as frozen_file:
    LEGACY_REFERENCE = {key: frozen_file[key] for key in frozen_file.files}


def test_strain_stress_heating_matches_legacy():
    """Strains match the classic results to machine precision, and stresses and heating do once the classic
    isotropic stress lambda tr(eps) is replaced by the kernel's (y2 - 2 mu dy1/dr) U.

    The two isotropic forms agree for a degree-2 harmonic (test_3d_fixes_01 checks it); the frozen potential rows
    are random, so here the classic stresses are converted before comparing.
    """
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
        s_ref = np.array(LEGACY_REFERENCE["classic_stresses"][sample_i], dtype=np.complex128)
        lame = bulk - (2.0 / 3.0) * shear
        # y2 = lame X + 2 mu dy1/dr with X = dy1/dr + (2 y1 - l(l+1) y3) / r, at l = 2.
        dy1_dr = (y[1] - lame * (2.0 * y[0] - 6.0 * y[2]) / radius) / (lame + 2.0 * shear)
        s_ref[:3] += (y[1] - 2.0 * shear * dy1_dr) * pot6[0] - lame * e_ref[:3].sum()
        h_ref = volumetric_heating(s_ref, np.ascontiguousarray(e_ref, dtype=np.complex128), frequency)
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
        # Against the sample's largest stress: the isotropic term is a difference of much larger terms, so a small
        # normal component carries their round-off.
        worst_s = max(worst_s, float(np.max(np.abs(s_c - s_ref)) / (np.max(np.abs(s_ref)) + 1e-30)))
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
