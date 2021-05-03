from .liquid_liquid import both_static as interface_LSt_LSt, both_dynamic as interface_LDy_LDy
from .liquid_solid import both_static as interface_LSt_SSt, both_dynamic as interface_LDy_SDy, \
    static_dynamic as interface_LSt_SDy, dynamic_static as interface_LDy_SSt
from .solid_liquid import both_static as interface_SSt_LSt, both_dynamic as interface_SDy_LDy, \
    static_dynamic as interface_SSt_LDy, dynamic_static as interface_SDy_LSt
from .solid_solid import both_static as interface_SSt_SSt, both_dynamic as interface_SDy_SDy, \
    static_dynamic as interface_SSt_SDy, dynamic_static as interface_SDy_SSt
