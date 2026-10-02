"""Materials the tests build again and again, shared so each test file states only its numbers.

Importable from any test (``from shared_materials import constant_solid``): pytest puts this folder on the import path
for Tests/conftest.py.
"""
from TidalPy.Material import Material, Phase


def constant_solid(
        density,
        bulk_modulus=None,
        shear_modulus=None,
        shear_viscosity=None,
        bulk_viscosity=None,
        eos_extra=None,
        **phase_kwargs):
    """A solid-only material whose every law is constant.

    Parameters
    ----------
    density : float
        Reference density [kg m-3].
    bulk_modulus, shear_modulus : float, optional
        Moduli [Pa]; a law left out is left out of the phase.
    shear_viscosity, bulk_viscosity : float, optional
        Viscosities [Pa s]; a law left out is left out of the phase.
    eos_extra : dict, optional
        More keys for the equation of state (thermal_expansion_1_k, reference_temperature_k, ...).
    **phase_kwargs
        The phase's other arguments (shear_rheology, thermal_conductivity_w_mk, ...).
    """
    eos = {"model": "constant", "reference_density_kg_m3": density}
    if bulk_modulus is not None:
        eos["bulk_modulus_pa"] = bulk_modulus
    eos.update(eos_extra or {})
    laws = {}
    if shear_modulus is not None:
        laws["shear_modulus"] = {"model": "constant", "shear_modulus_pa": shear_modulus}
    if shear_viscosity is not None:
        laws["shear_viscosity"] = {"model": "constant", "reference_viscosity_pas": shear_viscosity}
    if bulk_viscosity is not None:
        laws["bulk_viscosity"] = {"model": "constant", "reference_viscosity_pas": bulk_viscosity}
    return Material(solid=Phase(eos=eos, **laws, **phase_kwargs))
