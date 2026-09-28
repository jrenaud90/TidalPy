"""Tests that every demo and benchmark notebook runs top to bottom without raising (outputs are not compared).

Needs ``nbclient``, ``nbformat``, and ``ipykernel`` (the ``docs`` extra); without them the tests are skipped.
"""
from pathlib import Path

import pytest

nbformat = pytest.importorskip("nbformat")
nbclient = pytest.importorskip("nbclient")
pytest.importorskip("ipykernel")

REPO_ROOT = Path(__file__).resolve().parents[2]
NOTEBOOKS = sorted(
    path for path in (
        list((REPO_ROOT / "Demos").rglob("*.ipynb")) +
        list((REPO_ROOT / "Benchmarks" / "RadialSolver").glob("*.ipynb")) +
        list((REPO_ROOT / "Benchmarks" / "EOS").glob("*.ipynb")) +
        list((REPO_ROOT / "Benchmarks" / "Tides").glob("*.ipynb")))
    # Jupyter's autosave copies are not part of the repository.
    if ".ipynb_checkpoints" not in path.parts)

# Packages a notebook imports that TidalPy does not require; a notebook that needs a missing one is skipped.
OPTIONAL_PACKAGES = {"EOS_vs_BurnMan.ipynb": ("burnman",)}

# Generous: the slowest notebooks integrate an orbit or fill 3D grids, and an overloaded machine runs slower.
CELL_TIMEOUT = 900  # [s]

SETUP_TEMPLATE = """\
import os
os.environ["TIDALPY_TEST_MODE"] = "1"
from TidalPy.Structures.configs import worldpack
worldpack.get_worlds_dir = lambda: {worlds_dir!r}
"""


@pytest.mark.parametrize("notebook_path", NOTEBOOKS, ids=[path.name for path in NOTEBOOKS])
def test_notebook_executes(notebook_path, tmp_path):
    """The notebook executes in a fresh kernel from its own folder."""
    for package in OPTIONAL_PACKAGES.get(notebook_path.name, ()):
        pytest.importorskip(package)

    notebook = nbformat.read(notebook_path, as_version=4)
    worlds_dir = tmp_path / "Worlds"
    worlds_dir.mkdir()
    # The kernel is a separate process that Tests/conftest.py cannot reach, so an in-memory setup cell isolates it.
    notebook.cells.insert(0, nbformat.v4.new_code_cell(SETUP_TEMPLATE.format(worlds_dir=str(worlds_dir))))

    client = nbclient.NotebookClient(
        notebook,
        timeout=CELL_TIMEOUT,
        kernel_name="python3",
        resources={"metadata": {"path": str(notebook_path.parent)}})
    client.execute()
