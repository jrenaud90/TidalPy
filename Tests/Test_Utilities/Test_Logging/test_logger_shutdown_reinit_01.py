"""Tests that shutdown_logger silences the C++ logger and that the next init_logger reaches every extension again.

spdlog writes to the console directly, out of reach of ``caplog`` and ``capfd``, so the file sink is observed instead.
"""
import os
import subprocess
import sys

import pytest

from TidalPy.initialize import build_logging_config
from TidalPy.Utilities.logging import flush_logger, init_logger, log_warning, set_log_level, shutdown_logger


@pytest.fixture
def log_path(tmp_path):
    """A log file for the test; the package logging configuration is restored afterward."""
    yield tmp_path / "tidalpy_shutdown.log"
    init_logger(build_logging_config())


def configure_file_sink(path):
    init_logger({
        "console_level": "off",
        "file_level":    "trace",
        "log_to_file":   True,
        "log_file_path": str(path),
    })


def read_log(path):
    flush_logger()
    return path.read_text(encoding="utf-8") if path.exists() else ""


@pytest.mark.parametrize("turn_back_on", ["init_logger", "set_log_level"])
def test_shutdown_silences_until_turned_back_on(log_path, turn_back_on):
    """Nothing is written between shutdown_logger and the call that turns the logger back on."""
    configure_file_sink(log_path)
    shutdown_logger()
    log_warning("written while shut down")
    if turn_back_on == "init_logger":
        configure_file_sink(log_path)
    else:
        set_log_level("trace")
    log_warning("written after turning back on")
    text = read_log(log_path)
    assert "written while shut down" not in text
    assert "written after turning back on" in text


# The rheology extension logs the schema-patch line from C++ through its own copy of the logger pointer, which it
# wires when it is imported: here while the logger is shut down.
REINIT_SCRIPT = """
import sys
import TidalPy
from TidalPy.Utilities.logging import flush_logger, init_logger, shutdown_logger
already_imported = "TidalPy.Rheology.rheology" in sys.modules
shutdown_logger()
from TidalPy.Rheology import make_rheology
init_logger({
    "console_level": "off",
    "file_level": "trace",
    "log_to_file": True,
    "log_file_path": sys.argv[1],
})
model_path = sys.argv[2]
make_rheology("maxwell").save_binary(model_path)
with open(model_path, "rb") as model_file:
    data = bytearray(model_file.read())
data[6] = (data[6] + 1) % 256
with open(model_path, "wb") as model_file:
    model_file.write(bytes(data))
make_rheology("maxwell").load_binary(model_path)
flush_logger()
print("already imported" if already_imported else "imported while shut down")
"""


def test_extension_imported_while_shut_down_logs_after_reinit(tmp_path):
    """An extension first imported after shutdown_logger still logs once init_logger turns the logger back on."""
    subprocess_log_path = tmp_path / "tidalpy_reinit.log"
    environment = dict(os.environ)
    environment["TIDALPY_TEST_MODE"] = "1"
    result = subprocess.run(
        [sys.executable, "-c", REINIT_SCRIPT, str(subprocess_log_path), str(tmp_path / "maxwell.tpyb")],
        capture_output=True,
        text=True,
        timeout=300,
        cwd=tmp_path,
        env=environment,
    )
    assert result.returncode == 0, f"subprocess failed:\n{result.stderr}"
    if "already imported" in result.stdout:
        pytest.skip("importing TidalPy already imports the rheology extension, so the case cannot be set up")
    text = subprocess_log_path.read_text(encoding="utf-8") if subprocess_log_path.exists() else ""
    assert "schema patch differs" in text
