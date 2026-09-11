"""Old-versus-new radial solver comparison tests.

NOTE (0.9.0): every test here compares the `_x` solver against the classic TidalPy.RadialSolver and
loses its reference when the legacy tree is removed. Freeze the classic results as data (as
test_compare_benchmark_targets.py already does) or drop the test at that point.
"""
