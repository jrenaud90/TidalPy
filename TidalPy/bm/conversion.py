""" Common conversions between TidalPy and BurnMan """

burnman_property_name_conversion = {
    'thermal_expansion': 'thermal_expansivity',
    'grueneisen': 'grueneisen_parameter',
    'bulk_modulus': 'adiabatic_bulk_modulus',
    'heat_capacity': 'molar_heat_capacity_p',
}

burnman_property_value_conversion = {
    # Need to change BurnMan's molar heat capacity to specific
    'heat_capacity': 'molar'
}