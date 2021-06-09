#
# def build_diffeq(object_1, object_2):
#
#     # Check if object is planet or a dict. Pull out info based on that.
#
#     num_layers_obj1 = len(object_1)
#     thermal_conductivity_obj1 = tuple([layer.thermal_conductivity for layer in object_1])
#     thermal_diffusivity_obj1 = tuple([layer.thermal_diffusivity for layer in object_1])
#     static_shear_obj1 = tuple([layer.static_shear_modulus for layer in object_1])
#     static_viscosity_obj1 = tuple([layer.static_viscosity for layer in object_1])
#     partial_melt_func_obj1 = tuple([layer.rheology.partial_melting_model.func for layer in object_1])
#     viscosity_func_obj1 = tuple([layer.rheology.viscosity_model.func for layer in object_1])
#
#     @njit()
#     def diffeq(time, variables):
#
#         # TODO: orbital variation wil lbe first few variables
#         temperature_obj1 = variables[0:num_layers_obj1+1]
#         temperature_obj2 = variables[num_layers_obj1 + 1:num_layers_obj2 + 1]
#
#         # Update temperature-related items
#         viscosity_obj1 = list()
#         viscosity_obj2 = list()
#         for object_i, num_layers in enumerate((num_layers_obj1, num_layers_obj2)):
#             for layer_i in range(0, num_layers):
#                 pre_melt_viscosity =
#                 pre_melt_shear_modulus = static_shear[object_1][layer_i]