"""Benchmark task modules.

Importing this package imports every task module below, which registers its tasks with the
harness through the ``@benchmark`` decorator. ``run_all.py`` imports this package and then calls
``harness.run_suite()``. The task set mirrors the common operations taught in ``Demos_x/`` so the
tracked timings reflect real usage.
"""
from . import common_tasks  # noqa: F401
