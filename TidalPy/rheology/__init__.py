# Rheology Imports
from .rheology import Rheology as Rheology

from TidalPy.rheology.models import find_rheology as find_rheology
from TidalPy.rheology.models import Elastic as Elastic
from TidalPy.rheology.models import Newton as Newton
from TidalPy.rheology.models import Maxwell as Maxwell
from TidalPy.rheology.models import Voigt as Voigt
from TidalPy.rheology.models import Burgers as Burgers
from TidalPy.rheology.models import Andrade as Andrade
from TidalPy.rheology.models import SundbergCooper as SundbergCooper

# Alias rheologies
Sundberg = SundbergCooper
