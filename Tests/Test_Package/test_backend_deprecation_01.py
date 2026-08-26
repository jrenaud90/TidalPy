"""The package announces the classic-to-`_x` backend transition once per session at import time.

The warning uses the named category ``TidalPy.exceptions.TidalPyDeprecationWarning`` (a ``FutureWarning``
subclass, so it is visible by default) and points at the porting guide. A subprocess is used because the
pytest process has already imported TidalPy; the subprocess runs from a neutral working directory so the
installed package is imported rather than the repo source tree.
"""
import subprocess
import sys


def test_import_emits_backend_deprecation_warning(tmp_path):
    """Importing TidalPy emits exactly one TidalPyDeprecationWarning pointing at the porting guide."""
    code = (
        "import warnings\n"
        "with warnings.catch_warnings(record=True) as caught:\n"
        "    warnings.simplefilter('always')\n"
        "    import TidalPy\n"
        "from TidalPy.exceptions import TidalPyDeprecationWarning\n"
        "assert issubclass(TidalPyDeprecationWarning, FutureWarning)\n"
        "hits = [w for w in caught if issubclass(w.category, TidalPyDeprecationWarning)]\n"
        "assert len(hits) == 1, f'expected exactly one backend warning, got {len(hits)}'\n"
        "message = str(hits[0].message)\n"
        "assert 'future_structure' in message\n"
        "assert 'filterwarnings' in message\n"
        "print('BACKEND_WARNING_OK')\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        timeout=300,
        cwd=tmp_path,
    )
    assert result.returncode == 0, f"subprocess failed:\n{result.stderr}"
    assert "BACKEND_WARNING_OK" in result.stdout
