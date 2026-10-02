"""NumPy functions the tests use whose name depends on the NumPy version.

Importable from any test (``from numpy_compat import trapezoid``): pytest puts this folder on the import path for
Tests/conftest.py.
"""
import numpy as np

# The trapezoid rule: ``np.trapezoid`` from NumPy 2.0, ``np.trapz`` before it (TidalPy supports NumPy from 1.22).
trapezoid = getattr(np, "trapezoid", None) or np.trapz
