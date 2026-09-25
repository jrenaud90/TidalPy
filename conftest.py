"""Pytest bootstrap: import the installed (compiled) TidalPy, not the source tree.

The repository root contains the ``TidalPy`` source package alongside in-place build artifacts that can
be stale relative to the wheel installed into the active environment (``uv pip install .``). Pytest would
otherwise prepend the repo root to ``sys.path`` and import that source copy, so the ``_x`` extension
modules (and freshly added symbols such as ``update_constants_x``) would be missing or out of date. Drop
the repo root here so ``import TidalPy`` resolves to site-packages.
"""
import os
import sys

_REPO_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path[:] = [p for p in sys.path if os.path.abspath(p or os.getcwd()) != _REPO_ROOT]
