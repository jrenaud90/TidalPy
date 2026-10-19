"""Tests that every demo and benchmark notebook runs top to bottom without raising (outputs are not compared).

Needs ``nbclient``, ``nbformat``, and ``ipykernel`` (the ``docs`` extra); without them the tests are skipped.
"""
import pytest

from demo_notebooks import NOTEBOOKS

nbformat = pytest.importorskip("nbformat")
nbclient = pytest.importorskip("nbclient")
pytest.importorskip("ipykernel")

# Packages a notebook imports that TidalPy does not require; a notebook that needs a missing one is skipped.
OPTIONAL_PACKAGES = {"EOS_vs_BurnMan.ipynb": ("burnman",)}

# Generous: the slowest notebooks integrate an orbit or fill 3D grids, and an overloaded machine runs slower.
CELL_TIMEOUT = 900  # [s]

# Notebooks that cannot finish inside CELL_TIMEOUT, each with the reason, so a plain `pytest Tests/` passes without a
# command-line filter and the skip shows in the report.
TOO_SLOW = {
    "S02_thermal_orbital_evolution.ipynb": "integrates 5 Gyr of a coupled thermal-orbital system with System.evolve, "
                                           "about three minutes, longer than the rest of the suite together",
}

SETUP_TEMPLATE = """\
import os
os.environ["TIDALPY_DATA_DIR"] = {data_dir!r}
os.environ["TIDALPY_TEST_MODE"] = "1"
"""


# Notebooks that take more than about four seconds in the suite; they run only in the full suite.
SLOW = {
    "RadialSolver - General Benchmarks and Performance.ipynb", "P15_dynamic_liquid_stability.ipynb",
    "Homogeneous_Viscoelastic_Love_Numbers.ipynb", "Renaud2021_Dual_Body_Eccentric.ipynb",
    "Enceladus_Tobie_Roberts.ipynb", "P11_tidal_truncations.ipynb", "P08_tidal_maps_3d.ipynb",
    "Earth_Love_Numbers.ipynb",
}


@pytest.mark.parametrize(
    "notebook_path",
    [pytest.param(path, marks=pytest.mark.slow) if path.name in SLOW else path for path in NOTEBOOKS],
    ids=[path.name for path in NOTEBOOKS])
def test_notebook_executes(notebook_path, tmp_path):
    """The notebook executes in a fresh kernel from its own folder."""
    if notebook_path.name in TOO_SLOW:
        pytest.skip(f"{notebook_path.name} {TOO_SLOW[notebook_path.name]}")
    for package in OPTIONAL_PACKAGES.get(notebook_path.name, ()):
        pytest.importorskip(package)

    notebook = nbformat.read(notebook_path, as_version=4)
    # The kernel is a separate process, so a first cell gives it its own fresh data directory before it imports TidalPy
    # (as Tests/conftest.py does for the test process).
    setup = SETUP_TEMPLATE.format(data_dir=str(tmp_path / "data"))
    notebook.cells.insert(0, nbformat.v4.new_code_cell(setup))

    client = nbclient.NotebookClient(
        notebook,
        timeout=CELL_TIMEOUT,
        kernel_name="python3",
        resources={"metadata": {"path": str(notebook_path.parent)}})
    client.execute()
