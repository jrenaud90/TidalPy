"""The demo and benchmark notebooks the Test_Demos tests check (Jupyter's autosave copies left out).

Importable from any test in this folder (``from demo_notebooks import NOTEBOOKS``): pytest puts the folder on the
import path.
"""
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
NOTEBOOKS = sorted(
    path for path in (
        list((REPO_ROOT / "Demos").rglob("*.ipynb")) +
        list((REPO_ROOT / "Benchmarks" / "RadialSolver").glob("*.ipynb")) +
        list((REPO_ROOT / "Benchmarks" / "EOS").glob("*.ipynb")) +
        list((REPO_ROOT / "Benchmarks" / "Tides").glob("*.ipynb")))
    if ".ipynb_checkpoints" not in path.parts)
