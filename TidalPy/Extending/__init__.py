extended_configs = dict()

from TidalPy.utilities.dictionary_utils import nested_merge
from TidalPy.Extending.burnman import burnman_installed

# Add any extension configs to the extended config dict which will be added to the tidalpy configs
if burnman_installed:
    from TidalPy.Extending.burnman.burnman_defaultc import default_burnman_configs
    extended_configs = nested_merge(extended_configs, default_burnman_configs)