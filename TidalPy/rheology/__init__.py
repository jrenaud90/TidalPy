# Rheology Imports
from .rheology import Rheology
from TidalPy.rheology.models import (
    Elastic,
    Newton,
    Maxwell,
    Voigt,
    Burgers,
    Andrade,
    SundbergCooper)
from TidalPy.rheology.models import find_rheology

# Alias rheologies
Sundberg = SundbergCooper
