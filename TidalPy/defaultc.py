# These default configurations are saved to the user"s Application Data directory when TidalPy is first called.
# After the initial call, TidalPy will no longer look at this file.
# Config changes should be made in the AppData toml file.

default_config_str = """
[pathing]
    # Determine TidalPy"s output directory structure and naming scheme.
    save_directory = "TidalPy-Run"
    append_datetime = true

[debug]
    # Additional logging is used to track down bugs or issues. There could be a performance penalty in using this.
    extensive_logging = true
    # Additional numerical sanity checks will be performed. There could be a performance penalty in using this.
    extensive_checks = false

[logging]
    # Are TidalPy logs are stored to the current working directory or the default TidalPy data path.
    use_cwd = true

    # Should TidalPy write its debug/info log to `save_directory`?
    write_log_to_disk = false

    # Logging levels
    # Options: DEBUG, INFO, WARNING, ERROR
    file_level = "DEBUG"
    console_level = "INFO"
    console_error_level = "WARNING"

    # Determine if log should be printed to consol or written to disk when run in a jupyter notebook.
    print_log_notebook = false
    write_log_notebook = false

[configs]
    # Save these configurations to local run directory.
    save_configs_locally = false
    use_cwd_for_config = false
    use_cwd_for_world_dir = false

    # Determine how locally saved configs are handles
    overwrite_configs = true

[numba]
    # TODO: Remove this section if/when numba support is deprecated.
    # numba.njit can speed up many functions, but it also makes debugging and error tracing more difficult.
    #  If you are having problems try setting this to false.
    use_numba = true
    cache_numba = true

[tides]
    [tides.modes]
        # Set the minimum forcing frequency that is treated as zero: `|w| <= minimum_frequency --> w = 0`
        # The default value corresponds to a forcing period of around a billion years.
        minimum_frequency = 1.0e-17
    
    [tides.models]
        [tides.models.base]
            eccentricity_truncation_lvl = 6
            max_tidal_order_l = 2
            obliquity_tides_on = true
            use_planet_params_for_love_calc = true
            multiply_modes_by_sign = true
            slices = 100
        [tides.models.global_approx]
            eccentricity_truncation_lvl = 6
            max_tidal_order_l = 2
            obliquity_tides_on = true
            use_planet_params_for_love_calc = true
            multiply_modes_by_sign = true
            slices = 100
            fixed_q = 1000.0
            static_k2 = 0.37
            use_ctl = false
            # ctl_calc_method and fixed_dt used for CTL method
            ctl_calc_method = "linear_simple"
            fixed_dt = -7.27220521664304e-08  # (1. / 100.) * (2. * np.pi / (86400. * 10.))**(-1)
        [tides.models.layered]
            eccentricity_truncation_lvl = 6
            max_tidal_order_l = 2
            obliquity_tides_on = true
            use_planet_params_for_love_calc = false
            multiply_modes_by_sign = true
        [tides.models.multilayer]
            eccentricity_truncation_lvl = 6
            max_tidal_order_l = 2
            obliquity_tides_on = true
            multiply_modes_by_sign = true

[layers]

    # Known layer types
    [layers.ice]
        # General Layer Information
        type = "ice"
        # Switches
        use_surface_gravity = true
        use_bulk_density = true
        use_tidal_vol_frac = true
        is_tidally_active = true
        use_pressure_in_strength_calc = false

        # Burnman Information
        interp_temperature_range = [20.0, 273.15]

        # Material Information
        material = false
        material_source = false
        slices = 40

        # TODO: Eventually this material information should be offloaded to a material class via BM or BM Wrapper
        stefan = 0.0
        boundary_temperature_ratio = 1.0
        heat_fusion = 2.84e5
        shear_modulus = 9.2e9
        thermal_conductivity = 2.3
        thermal_diffusivity = 1.15e-6
        thermal_expansion = 5.0e-5

        [layers.ice.rheology]
            model = "off"
            voigt_compliance_offset = 0.2
            voigt_viscosity_offset = 0.02
            alpha = 0.3333333333333333  # 1/3
            zeta = 1.0
            critical_freq = 7.27220521664303e-07  # 2 * pi / (86400 * 100)
            critical_freq_falloff = 30.0

        [layers.ice.partial_melting]
            model = "off"
            solidus = 270.0
            liquidus = 273.15
            liquid_shear = 1.0e-5
        
        [layers.ice.liquid_viscosity]
            model = "reference"
            reference_viscosity = 0.89e-3  # [Pa s]
            reference_temperature = 298.15  # [K]
            molar_activation_energy = 1.62e4
            molar_activation_volume = 0.0
        
        [layers.ice.solid_viscosity]
            model = "arrhenius"
            # Matching Moore2006 for Volume Diffusion
            arrhenius_coeff = 1.1037527593819e07
            additional_temp_dependence = true
            stress = 1.0
            stress_expo = 1.0
            grain_size = 5.0e-4
            grain_size_expo = 2.0
            molar_activation_energy = 59.4e3
            # FIXME: I get crazy low viscosity values when I have this activation volume. For now not assuming pressure dependence.
            # molar_activation_volume = -1.3e-5
            molar_activation_volume = 0.0

        [layers.ice.radiogenics]
            model = "off"
            use_full_layer_mass_for_radiogenics = true
            radiogenic_layer_mass_fraction = 1.0
        
        [layers.ice.cooling]
            model = "convection"
            convection_alpha = 1.0
            convection_beta = 0.3333333333333333  # 1/3
            critical_rayleigh = 1600.0

    [layers.rock]
        # General Layer Information
        type = "rock"
        # Switches
        use_surface_gravity = true
        use_bulk_density = true
        use_tidal_vol_frac = true
        is_tidally_active = true
        use_pressure_in_strength_calc = false

        # Burnman Information
        interp_temperature_range = [100.0, 3200.0]

        # Material Information
        material = false
        material_source = false
        slices = 40

        # TODO: Eventually this material information should be offloaded to a material class via BM or BM Wrapper
        stefan = 0.0
        boundary_temperature_ratio = 1.0
        heat_fusion = false
        shear_modulus = 6.0e10
        thermal_conductivity = 3.75
        thermal_diffusivity = 9.15751e-7
        thermal_expansion = 5.2e-5

        [layers.rock.rheology]
            model = "maxwell"
            voigt_compliance_offset = 0.2
            voigt_viscosity_offset = 0.02
            alpha = 0.3333333333333333  # 1/3
            zeta = 1.0
            critical_freq = 7.27220521664303e-07  # 2 * pi / (86400 * 100)
            critical_freq_falloff = 30.0

        [layers.rock.partial_melting]
            model = "henning"
            solidus = 1600.0
            liquidus = 2000.00
            liquid_shear = 1.0e-5
            fs_visc_power_slope = 27000.0
            fs_visc_power_phase = 1.0
            fs_shear_power_slope = 82000.0
            fs_shear_power_phase = 40.6
            crit_melt_frac = 0.5
            crit_melt_frac_width = 0.05
            hn_visc_slope_1 = 13.5
            hn_visc_slope_2 = 370.0
            hn_shear_param_1 = 40000.0
            hn_shear_param_2 = 25.0
            hn_shear_falloff_slope = 700.0
        
        [layers.rock.liquid_viscosity]
            model = "reference"
            reference_viscosity = 0.2  # [Pa s]
            reference_temperature = 2000.0  # [K]
            molar_activation_energy = 6.64e-20
            molar_activation_volume = 0.0
        
        [layers.rock.solid_viscosity]
            model = "reference"
            reference_viscosity = 1.0e22
            reference_temperature = 1000.0
            molar_activation_energy = 300000.0
            molar_activation_volume = 0.0
        
        [layers.rock.radiogenics]
            model = "isotope"
            # isotopes can either be a string that references one of the pre-defined sets above,
            #  or a custom set (use the format found for [physics.radiogenics]).
            isotopes = "modern_day_chondritic"
            radiogenic_layer_mass_fraction = 1.0
        
        [layers.rock.cooling]
            model = "convection"
            convection_alpha = 1.0
            convection_beta = 0.3333333333333333  # 1/3
            critical_rayleigh = 1100.0

    [layers.iron]
        # General Layer Information
        type = "iron"
        # Switches
        use_surface_gravity = true
        use_bulk_density = true
        use_tidal_vol_frac = true
        is_tidally_active = false
        use_pressure_in_strength_calc = false

        # Burnman Information
        interp_temperature_range = [500.0, 5000.0]

        # Material Information
        material = false
        material_source = false
        slices = 40

        # TODO: Eventually this material information should be offloaded to a material class via BM or BM Wrapper
        stefan = 0.0
        boundary_temperature_ratio = 1.1
        heat_fusion = false
        shear_modulus = 5.25e10
        thermal_conductivity = 7.95
        thermal_diffusivity = 1.98949e-6
        thermal_expansion = 1.2e-5

        [layers.iron.rheology]
            model = "maxwell"
            voigt_compliance_offset = 0.2
            voigt_viscosity_offset = 0.02
            alpha = 0.3333333333333333  # 1/3
            zeta = 1.0
            critical_freq = 7.27220521664303e-07  # 2 * pi / (86400 * 100)
            critical_freq_falloff = 30.0
        
        [layers.iron.partial_melting]
            model = "off"
            solidus = 4000.0
            liquidus = 5000.00
            liquid_shear = 1.0e-5
        
        [layers.iron.liquid_viscosity]
            # These values match Wijs et al 1998 (their work actually does not show much change in the liquid visc
            #    at Earth"s core pressure, so a constant model may not be too incorrect).
            model = "constant"
            reference_viscosity = 1.3e-2  # [Pa s]

        [layers.iron.solid_viscosity]
            model = "constant"
            reference_viscosity = 1.0e20
        
        [layers.iron.radiogenics]
            model = "off"
            use_full_layer_mass_for_radiogenics = true
            radiogenic_layer_mass_fraction = 1.0

        [layers.iron.cooling]
            model = "off"
    
[worlds]
    # Should TidalPy save world information to `save_directory`?
    save_worlds_to_disk = true
    autosave_worlds_dir = true

    # Known world types
    [worlds.types]
        [worlds.types.base]
            name = "unknown_world_base_type"
            store_tides_config_in_world = true
            force_spin_sync = true
            equilibrium_insolation_model = "williams" # Options are no_eccentricity, williams, or mendez
            fraction_internal_heating_to_surface = 1.0
            emissivity = 0.9
            albedo = 0.3
            use_real_moi = true
            surface_pressure = 0.0
            slices = 40
        
        [worlds.types.simple_tidal]
            name = "unknown_world_simple_tide_type"
            store_tides_config_in_world = true
            force_spin_sync = true
            equilibrium_insolation_model = "williams" # Options are no_eccentricity, williams, or mendez
            fraction_internal_heating_to_surface = 1.0
            emissivity = 0.9
            albedo = 0.3
            use_real_moi = true
            tides_on = true
            surface_pressure = 0.0
            slices = 40
        
        [worlds.types.gas_giant]
            name = "unknown_world_gas_giant_type"
            store_tides_config_in_world = true
            force_spin_sync = true
            equilibrium_insolation_model = "williams" # Options are no_eccentricity, williams, or mendez
            fraction_internal_heating_to_surface = 1.0
            emissivity = 0.9
            albedo = 0.3
            use_real_moi = true
            tides_on = true
            surface_pressure = 0.0
            slices = 40
        
        [worlds.types.star]
            name = "unknown_world_star_type"
            store_tides_config_in_world = true
            force_spin_sync = true
            equilibrium_insolation_model = "williams" # Options are no_eccentricity, williams, or mendez
            fraction_internal_heating_to_surface = 1.0
            emissivity = 0.9
            albedo = 0.3
            use_real_moi = true
            tides_on = false
            surface_pressure = 0.0
            slices = 40
        
        [worlds.types.layered]
            name = "unknown_world_layered_type"
            store_tides_config_in_world = true
            force_spin_sync = true
            equilibrium_insolation_model = "williams" # Options are no_eccentricity, williams, or mendez
            fraction_internal_heating_to_surface = 1.0
            emissivity = 0.9
            albedo = 0.3
            use_real_moi = true
            tides_on = true
            surface_pressure = 0.0
            slices = 40

[physics]
    [physics.radiogenics]
        [physics.radiogenics.known_isotope_data]
            [physics.radiogenics.known_isotope_data.modern_day_chondritic]
                ref_time = 4600.0

                # Based off Hussmann & Spohn 2004 and Turcotte & Schubert 2001
                [physics.radiogenics.known_isotope_data.modern_day_chondritic.U238]
                    iso_mass_fraction = 0.9928
                    hpr = 9.48e-5
                    half_life = 4470.0
                    element_concentration = 0.012e-6
                [physics.radiogenics.known_isotope_data.modern_day_chondritic.U235]
                    iso_mass_fraction = 0.0071
                    hpr = 5.69e-4
                    half_life = 704.0
                    element_concentration = 0.012e-6
                [physics.radiogenics.known_isotope_data.modern_day_chondritic.Th232]
                    iso_mass_fraction = 0.9998
                    hpr = 2.69e-5
                    half_life = 14000.0
                    element_concentration = 0.04e-6
                [physics.radiogenics.known_isotope_data.modern_day_chondritic.K40]
                    iso_mass_fraction = 1.19e-4
                    hpr = 2.92e-5
                    half_life = 1250.0
                    element_concentration = 840.0e-6
            
            [physics.radiogenics.known_isotope_data.LLRI_and_SLRI]
                ref_time = 4600.0

                # Based off Castillo-Rogez et al 2007
                [physics.radiogenics.known_isotope_data.LLRI_and_SLRI.U238]
                    iso_mass_fraction = 0.9928
                    hpr = 9.465e-5
                    half_life = 4468.0
                    element_concentration = 0.026e-6
                [physics.radiogenics.known_isotope_data.LLRI_and_SLRI.U235]
                    iso_mass_fraction = 0.0071
                    hpr = 5.687e-4
                    half_life = 703.81
                    element_concentration = 0.0082e-6
                [physics.radiogenics.known_isotope_data.LLRI_and_SLRI.Th232]
                    iso_mass_fraction = 1.0
                    hpr = 2.638e-5
                    half_life = 14025.0
                    element_concentration = 0.0538e-6
                [physics.radiogenics.known_isotope_data.LLRI_and_SLRI.K40]
                    iso_mass_fraction = 1.176e-4
                    hpr = 2.917e-5
                    half_life = 1277.0
                    element_concentration = 1.104e-6
                [physics.radiogenics.known_isotope_data.LLRI_and_SLRI.Mn53]
                    iso_mass_fraction = 2.0e-5
                    hpr = 0.027
                    half_life = 3.7
                    element_concentration = 0.0257e-6
                [physics.radiogenics.known_isotope_data.LLRI_and_SLRI.Fe60]
                    iso_mass_fraction = 1.0e-6
                    hpr = 0.07
                    half_life = 1.5
                    element_concentration = 0.1e-6
                [physics.radiogenics.known_isotope_data.LLRI_and_SLRI.Al26]
                    iso_mass_fraction = 5.0e-5
                    hpr = 0.146
                    half_life = 0.72
                    element_concentration = 0.6e-6
"""
