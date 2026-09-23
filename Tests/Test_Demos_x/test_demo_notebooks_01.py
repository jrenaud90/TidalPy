"""Every new-backend demo and benchmark notebook runs top to bottom against the installed package.

Nothing else in the suite executes the notebooks, so an API change can break every demo while the tests stay green
(a layer constructor argument that moved onto the material EOS broke four notebooks that way). Each notebook is run
in a fresh kernel from its own folder. A setup cell, added to the in-memory copy only, puts TidalPy in test mode and
points the world pack at a temporary folder, as ``Tests/conftest.py`` does for this process: the kernel is a separate
process, so that fixture does not reach it. Outputs are not compared; a cell that raises fails the test.

Needs ``nbclient``, ``nbformat``, and ``ipykernel`` (the ``docs`` extra); without them the tests are skipped.
"""
from pathlib import Path

import pytest

nbformat = pytest.importorskip("nbformat")
nbclient = pytest.importorskip("nbclient")
pytest.importorskip("ipykernel")

REPO_ROOT = Path(__file__).resolve().parents[2]
NOTEBOOKS = sorted(
    list((REPO_ROOT / "Demos_x").rglob("*.ipynb")) +
    list((REPO_ROOT / "Benchmarks_x" / "RadialSolver").glob("*.ipynb")) +
    list((REPO_ROOT / "Benchmarks_x" / "EOS").glob("*.ipynb")))

# Packages a notebook imports that TidalPy does not require; a notebook that needs a missing one is skipped.
OPTIONAL_PACKAGES = {"EOS_vs_BurnMan.ipynb": ("burnman",)}

# Generous: the slowest notebooks integrate an orbit or fill 3D grids, and an overloaded machine runs slower.
CELL_TIMEOUT = 900  # [s]

SETUP_TEMPLATE = """\
import os
os.environ["TIDALPY_TEST_MODE"] = "1"
from TidalPy.structures_x.configs import worldpack
worldpack.get_worlds_x_dir = lambda: {worlds_dir!r}
"""


@pytest.mark.parametrize("notebook_path", NOTEBOOKS, ids=[path.name for path in NOTEBOOKS])
def test_notebook_executes(notebook_path, tmp_path):
    for package in OPTIONAL_PACKAGES.get(notebook_path.name, ()):
        pytest.importorskip(package)

    notebook = nbformat.read(notebook_path, as_version=4)
    worlds_dir = tmp_path / "Worlds_x"
    worlds_dir.mkdir()
    notebook.cells.insert(0, nbformat.v4.new_code_cell(SETUP_TEMPLATE.format(worlds_dir=str(worlds_dir))))

    client = nbclient.NotebookClient(
        notebook,
        timeout=CELL_TIMEOUT,
        kernel_name="python3",
        resources={"metadata": {"path": str(notebook_path.parent)}})
    client.execute()
