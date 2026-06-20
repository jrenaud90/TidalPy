"""Default configuration for the new ``_x`` structure of TidalPy.

This module holds the default contents of ``TidalPy_Configs_x.toml``, the
configuration file for the rebuilt ``_x`` class system (schema ``0.2.0``). It is
written to the user's TidalPy data directory (next to the legacy
``TidalPy_Configs.toml``) on first use and is then user-editable.

The file has two roles:

1. ``[numerical]`` carries the numerical settings consumed by the C++ config
   singleton for ``_x`` code (see ``TidalPy.constants.update_constants_x``).
2. ``[layers.<type>]`` carries per-material default parameters keyed by a layer's
   material ``type`` (``gas``, ``mantle_rock``, ``ice``, ``hp_ice``, ``iron``).
   These supply the second tier of the world builder's default-resolution chain:
   a value the user omits is taken from the matching ``[layers.<type>]`` section
   here, and only if that is also absent does the C++/Cython constructor default
   apply.

IMPORTANT (for future work): any new default configuration that belongs to the
``_x`` system should be added here (and so to ``TidalPy_Configs_x.toml``), NOT to
the legacy ``defaultc.py`` / ``TidalPy_Configs.toml``.
"""

from TidalPy import version

# Schema version for the _x configuration / world format. Kept in sync with
# TidalPy.structures_x.configs.toml_loader.SCHEMA_VERSION.
SCHEMA_VERSION_X = "0.2.0"

default_config_x_str = f"""
schema_version = "{SCHEMA_VERSION_X}"

# =====================================================================================================================
# Numerical settings (consumed by the C++ config singleton for _x code)
# =====================================================================================================================
[numerical]
    # Forcing-frequency handling: |w| <= minimum_frequency is treated as zero.
    minimum_frequency = 1.0e-14
    maximum_frequency = 1.0e8
    min_spin_orbit_diff = 1.0e-10
    # Material floors.
    minimum_viscosity = 100.0
    minimum_modulus = 1.0e-3
    # Geometry floor: a layer thinner than this is ignored.
    minimum_layer_thickness = 0.1
    # Debug helper.
    test_constant = 42.0


# =====================================================================================================================
# Global (1D) tidal-dissipation defaults
#
# Used by the world builder when a world's `[tides]` table omits a value. The default
# dissipation model is chosen per world family; the truncation/degree settings apply to
# the global tidal-mode solve. The per-degree analytic parameters (fixed_k/fixed_q/fixed_dt,
# lists indexed from degree l = 2) are only consumed by the analytic models (cpl/ctl/ctl_q)
# and are easily overridden per world.
# =====================================================================================================================
[tides]
    min_degree_l = 2
    max_degree_l = 2
    eccentricity_trunc_lvl = 6
    obliquity_trunc_lvl = "off"

    # Per-degree static potential Love numbers k_l (index 0 -> l=2). Falls off roughly as
    # k_l ~ k_2 / (l - 1) for a soft, near-homogeneous body (k_2 ~ 0.3).
    fixed_k = [0.3, 0.15, 0.1, 0.075, 0.06, 0.05, 0.0429, 0.0375, 0.0333]
    # Per-degree tidal quality factors Q_l (a single constant-Q value across degrees).
    fixed_q = [100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0]
    # Per-degree tidal time lags dt_l [s] (a single constant time lag across degrees;
    # ~600 s is an Earth-like value).
    fixed_dt = [600.0, 600.0, 600.0, 600.0, 600.0, 600.0, 600.0, 600.0, 600.0]

    # Default global dissipation model per world family.
    [tides.default_model]
        star = "fixed_q"
        gasgiant = "fixed_dt"
        terrestrial = "rheology"
        layered = "rheology"


# =====================================================================================================================
# Per-material layer defaults (keyed by a layer's material `type`)
#
# A layer in a world TOML names a `class` (base | physics | solidliquid | gas) and,
# optionally, a material `type` (one of the sections below). Any layer parameter or
# physics-model sub-table the user omits falls back to the matching section here.
# Model sub-tables / scalar keys that a given layer `class` cannot hold are ignored
# for that layer, so the same material defaults can be reused across layer classes.
# =====================================================================================================================

# Iron (metallic core). Typically a non-tidal solidliquid layer.
[layers.iron]
    shear_modulus_static_pa = 5.25e10
    bulk_modulus_static_pa = 1.6e11
    thermal_conductivity_ref_w_mk = 7.95
    thermal_expansion_ref_1_k = 1.2e-5
    heat_capacity_ref_j_kgk = 840.0
    solidus_temperature_k = 4000.0
    liquidus_temperature_k = 5000.0
    reference_density_kg_m3 = 8000.0
    reference_temperature_k = 4000.0

    [layers.iron.eos]
        model = "constant"
        reference_density_kg_m3 = 8000.0

    [layers.iron.shear_rheology]
        model = "maxwell"

    [layers.iron.bulk_rheology]
        model = "elastic"

    [layers.iron.shear_viscosity]
        model = "constant"
        reference_viscosity = 1.0e20

    [layers.iron.bulk_viscosity]
        model = "constant"
        reference_viscosity = 1.0e22

    [layers.iron.partial_melt]
        model = "off"
        solidus_k = 4000.0
        liquidus_k = 5000.0

    [layers.iron.cooling]
        model = "off"

    [layers.iron.radiogenics]
        model = "off"

# Silicate mantle rock. The canonical tidally active solid layer.
[layers.mantle_rock]
    shear_modulus_static_pa = 6.0e10
    bulk_modulus_static_pa = 2.0e11
    thermal_conductivity_ref_w_mk = 3.75
    thermal_expansion_ref_1_k = 5.2e-5
    heat_capacity_ref_j_kgk = 1200.0
    activation_energy_j_mol = 3.0e5
    activation_volume_m3_mol = 5.0e-6
    solidus_temperature_k = 1600.0
    liquidus_temperature_k = 2000.0
    reference_density_kg_m3 = 3500.0
    reference_temperature_k = 1600.0
    melt_viscosity_reduction = 25.0

    [layers.mantle_rock.eos]
        model = "constant"
        reference_density_kg_m3 = 3500.0

    [layers.mantle_rock.shear_rheology]
        model = "andrade"
        alpha = 0.3
        zeta = 1.0

    [layers.mantle_rock.bulk_rheology]
        model = "elastic"

    [layers.mantle_rock.shear_viscosity]
        model = "reference"
        reference_viscosity = 1.0e22
        reference_temperature = 1000.0
        molar_activation_energy = 3.0e5
        molar_activation_volume = 0.0

    [layers.mantle_rock.bulk_viscosity]
        model = "constant"
        reference_viscosity = 1.0e22

    [layers.mantle_rock.partial_melt]
        model = "henning"
        solidus_k = 1600.0
        liquidus_k = 2000.0
        liquid_shear_pa = 1.0e-5
        crit_melt_frac = 0.5
        crit_melt_frac_width = 0.05
        hn_visc_slope_1 = 13.5
        hn_visc_falloff_slope = 370.0
        hn_shear_param_1 = 40000.0
        hn_shear_param_2 = 25.0
        hn_shear_falloff_slope = 700.0

    [layers.mantle_rock.cooling]
        model = "convection"
        convection_alpha = 1.0
        convection_beta = 0.3333333333333333
        critical_rayleigh = 1100.0

    [layers.mantle_rock.radiogenics]
        model = "isotope"
        isotopes = "modern_day_chondritic"

# Low-pressure water ice (ice Ih). Tidally active outer-shell material.
[layers.ice]
    shear_modulus_static_pa = 3.3e9
    bulk_modulus_static_pa = 9.2e9
    thermal_conductivity_ref_w_mk = 2.3
    thermal_expansion_ref_1_k = 5.0e-5
    heat_capacity_ref_j_kgk = 2000.0
    solidus_temperature_k = 250.0
    liquidus_temperature_k = 273.15
    reference_density_kg_m3 = 1000.0
    reference_temperature_k = 250.0

    [layers.ice.eos]
        model = "constant"
        reference_density_kg_m3 = 1000.0

    [layers.ice.shear_rheology]
        model = "maxwell"

    [layers.ice.bulk_rheology]
        model = "elastic"

    [layers.ice.shear_viscosity]
        model = "arrhenius"
        arrhenius_coeff = 1.1037527593819e07
        additional_temp_dependence = true
        stress = 1.0
        stress_expo = 1.0
        grain_size = 5.0e-4
        grain_size_expo = 2.0
        molar_activation_energy = 59.4e3
        molar_activation_volume = 0.0

    [layers.ice.bulk_viscosity]
        model = "constant"
        reference_viscosity = 1.0e22

    [layers.ice.partial_melt]
        model = "off"
        solidus_k = 250.0
        liquidus_k = 273.15

    [layers.ice.cooling]
        model = "convection"
        convection_alpha = 1.0
        convection_beta = 0.3333333333333333
        critical_rayleigh = 1600.0

    [layers.ice.radiogenics]
        model = "off"

# High-pressure water ice (e.g. ice VI/VII in large icy worlds). Denser and stiffer.
[layers.hp_ice]
    shear_modulus_static_pa = 6.0e9
    bulk_modulus_static_pa = 1.4e10
    thermal_conductivity_ref_w_mk = 2.3
    thermal_expansion_ref_1_k = 4.0e-5
    heat_capacity_ref_j_kgk = 2000.0
    solidus_temperature_k = 270.0
    liquidus_temperature_k = 300.0
    reference_density_kg_m3 = 1300.0
    reference_temperature_k = 270.0

    [layers.hp_ice.eos]
        model = "constant"
        reference_density_kg_m3 = 1300.0

    [layers.hp_ice.shear_rheology]
        model = "andrade"
        alpha = 0.3
        zeta = 1.0

    [layers.hp_ice.bulk_rheology]
        model = "elastic"

    [layers.hp_ice.shear_viscosity]
        model = "arrhenius"
        arrhenius_coeff = 1.1037527593819e07
        additional_temp_dependence = true
        stress = 1.0
        stress_expo = 1.0
        grain_size = 5.0e-4
        grain_size_expo = 2.0
        molar_activation_energy = 59.4e3
        molar_activation_volume = 0.0

    [layers.hp_ice.bulk_viscosity]
        model = "constant"
        reference_viscosity = 1.0e22

    [layers.hp_ice.partial_melt]
        model = "off"
        solidus_k = 270.0
        liquidus_k = 300.0

    [layers.hp_ice.cooling]
        model = "convection"
        convection_alpha = 1.0
        convection_beta = 0.3333333333333333
        critical_rayleigh = 1600.0

    [layers.hp_ice.radiogenics]
        model = "off"

# Gaseous envelope (e.g. a gas-giant layer). Uses the `gas` layer class; the
# cooling/radiogenics sections below are ignored for that class.
[layers.gas]
    mean_molecular_weight_kg_mol = 2.22e-3
    adiabatic_index = 1.4
    reference_temperature_k = 165.0
    reference_density_kg_m3 = 1000.0
    shear_modulus_static_pa = 0.0
    bulk_modulus_static_pa = 1.0e5

    [layers.gas.eos]
        model = "constant"
        reference_density_kg_m3 = 1000.0

    [layers.gas.shear_rheology]
        model = "elastic"

    [layers.gas.bulk_rheology]
        model = "elastic"
"""
