from libcpp import bool as cpp_bool
        
cdef struct ODEInput:
    double frequency
    double G_to_use
    double degree_l
    void* eos_input_ptr
    int gravity_index
    int pressure_index
    cpp_bool update_bulk
    cpp_bool update_shear
