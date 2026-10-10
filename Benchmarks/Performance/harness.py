"""Lightweight performance-tracking harness for TidalPy.

Each run records a JSON file under ``results/`` that captures the wall-clock timing of a
set of registered tasks together with the environment they ran in: the date, the TidalPy
version, the git commit, the operating system, the CPU, and the Python version. The
companion ``Perf_Trends.ipynb`` notebook loads every ``results/*.json`` file and plots how
a task's timing has moved across those axes, so performance can be tracked as a function of
date, machine, and TidalPy version now and into the future.

The harness depends only on the standard library plus ``psutil`` so it runs anywhere
TidalPy itself runs. Register tasks with the ``@benchmark`` decorator, then call
``run_suite()`` (or run ``run_all.py``).
"""
from __future__ import annotations

import json
import platform
import statistics
import subprocess
import timeit
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Iterable, Optional

try:
    import psutil
except ImportError:
    psutil = None

RESULTS_DIR = Path(__file__).resolve().parent / "results"

# Version stamp written into every results JSON so future readers can detect layout changes.
RESULTS_SCHEMA_VERSION = "1.0"


# =====================================================================================================================
# Environment capture
# =====================================================================================================================
def _tidalpy_version() -> str:
    import TidalPy
    for attr in ("version", "__version__"):
        value = getattr(TidalPy, attr, None)
        if value:
            return str(value)
    try:
        from importlib.metadata import version
        return version("TidalPy")
    except Exception:
        return "unknown"


def _git_commit() -> str:
    repo_root = Path(__file__).resolve().parents[2]
    try:
        out = subprocess.run(
            ["git", "-C", str(repo_root), "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, timeout=10,
        )
        if out.returncode == 0:
            return out.stdout.strip()
    except Exception:
        pass
    return "unknown"


def capture_environment() -> dict:
    """Collect the machine, OS, version, and date metadata stamped onto every run."""
    cpu_max_mhz = None
    physical = logical = None
    if psutil is not None:
        try:
            freq = psutil.cpu_freq()
            if freq is not None:
                cpu_max_mhz = round(freq.max or freq.current, 1)
        except Exception:
            cpu_max_mhz = None
        physical = psutil.cpu_count(logical=False)
        logical = psutil.cpu_count(logical=True)
    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "tidalpy_version": _tidalpy_version(),
        "git_commit": _git_commit(),
        "os": platform.system(),
        "os_release": platform.release(),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "cpu_cores_physical": physical,
        "cpu_cores_logical": logical,
        "cpu_max_mhz": cpu_max_mhz,
        "python_version": platform.python_version(),
        "node": platform.node(),
    }


# =====================================================================================================================
# Task registry
# =====================================================================================================================
class Task:
    """A named benchmark: a zero-argument callable timed with a fixed ``timeit`` configuration.

    ``number`` is the inner ``timeit`` loop count; leave it ``None`` to auto-calibrate per task
    with ``timeit`` autoranging. ``repeats`` is how many independent timing samples are taken;
    the per-call minimum, median, mean, and standard deviation are recorded from those samples.
    ``setup`` runs once before timing (for example, to build a world the timed call then reuses).
    """

    def __init__(
        self,
        name: str,
        func: Callable[[], object],
        *,
        group: str = "general",
        repeats: int = 7,
        number: Optional[int] = None,
        setup: Optional[Callable[[], object]] = None,
        note: str = "",
    ):
        self.name = name
        self.func = func
        self.group = group
        self.repeats = repeats
        self.number = number
        self.setup = setup
        self.note = note


_REGISTRY: list[Task] = []


def benchmark(
    name: str,
    *,
    group: str = "general",
    repeats: int = 7,
    number: Optional[int] = None,
    setup: Optional[Callable[[], object]] = None,
    note: str = "",
) -> Callable[[Callable[[], object]], Callable[[], object]]:
    """Register a zero-argument callable as a benchmark task and return it unchanged."""
    def wrap(func: Callable[[], object]) -> Callable[[], object]:
        _REGISTRY.append(Task(name, func, group=group, repeats=repeats, number=number, setup=setup, note=note))
        return func
    return wrap


def registered_tasks() -> list[Task]:
    return list(_REGISTRY)


def clear_registry() -> None:
    _REGISTRY.clear()


# =====================================================================================================================
# Timing
# =====================================================================================================================
def _autorange(func: Callable[[], object]) -> int:
    number, _ = timeit.Timer(func).autorange()
    return number


def time_task(task: Task) -> dict:
    """Time one task and return per-call min/median/mean/stdev in seconds."""
    if task.setup is not None:
        task.setup()
    number = task.number or _autorange(task.func)
    raw = timeit.Timer(task.func).repeat(repeat=task.repeats, number=number)
    per_call = [sample / number for sample in raw]
    return {
        "name": task.name,
        "group": task.group,
        "note": task.note,
        "number": number,
        "repeats": task.repeats,
        "min_s": min(per_call),
        "median_s": statistics.median(per_call),
        "mean_s": statistics.fmean(per_call),
        "stdev_s": statistics.pstdev(per_call) if len(per_call) > 1 else 0.0,
    }


def _sanitize(name: str) -> str:
    keep = "-_."
    return "".join(ch if (ch.isalnum() or ch in keep) else "_" for ch in name)


def run_suite(
    tasks: Optional[Iterable[Task]] = None,
    *,
    output_dir: Path = RESULTS_DIR,
    label: str = "",
    verbose: bool = True,
) -> Path:
    """Run every task, then write one JSON record for the whole run and return its path."""
    task_list = list(tasks) if tasks is not None else registered_tasks()
    if not task_list:
        raise RuntimeError("No benchmark tasks registered. Import a task module before running the suite.")
    env = capture_environment()
    results = []
    for task in task_list:
        if verbose:
            print(f"  timing {task.name} ...", flush=True)
        # A failing task records an error row instead of discarding the whole run.
        try:
            results.append(time_task(task))
        except Exception as error:
            results.append({"name": task.name, "group": task.group, "error": repr(error)})
            if verbose:
                print(f"    task {task.name} failed: {error!r}", flush=True)
    payload = {"schema_version": RESULTS_SCHEMA_VERSION, "environment": env, "label": label, "results": results}

    output_dir.mkdir(parents=True, exist_ok=True)
    stamp = env["timestamp_utc"].replace(":", "").replace("-", "")
    label_part = f"_{_sanitize(label)}" if label else ""
    base_name = _sanitize(f"{stamp}_{env['node']}_{env['tidalpy_version']}") + label_part
    path = output_dir / f"{base_name}.json"
    # Never silently overwrite a same-second run; append a counter instead.
    counter = 1
    while path.exists():
        path = output_dir / f"{base_name}_{counter:02d}.json"
        counter += 1
    path.write_text(json.dumps(payload, indent=2))
    if verbose:
        print(f"Wrote {len(results)} task result(s) to {path}")
    return path


# =====================================================================================================================
# Result loading (used by the trend notebook)
# =====================================================================================================================
def load_results(output_dir: Path = RESULTS_DIR) -> list[dict]:
    """Flatten every results JSON into one row per (run, task) for trend analysis.

    Each row merges the run environment with that task's timing, so it can be dropped straight
    into a pandas DataFrame and grouped by ``name``, ``tidalpy_version``, ``os``, or date.
    """
    rows = []
    for path in sorted(output_dir.glob("*.json")):
        # A single malformed file should not take down the whole trend view.
        try:
            payload = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError) as error:
            print(f"Skipping unreadable results file {path.name}: {error!r}")
            continue
        env = payload.get("environment", {})
        label = payload.get("label", "")
        for res in payload.get("results", []):
            if "error" in res:
                continue
            rows.append({**env, **res, "run_label": label, "source_file": path.name})
    return rows


def load_results_df(output_dir: Path = RESULTS_DIR):
    """Return :func:`load_results` as a pandas DataFrame, one row per (run, task).

    ``pandas`` is imported lazily so the rest of the harness stays importable without it.
    """
    import pandas as pd

    frame = pd.DataFrame(load_results(output_dir))
    if not frame.empty:
        frame["timestamp_utc"] = pd.to_datetime(frame["timestamp_utc"])
        frame["min_ms"] = frame["min_s"] * 1.0e3
        frame["median_ms"] = frame["median_s"] * 1.0e3
        frame["stdev_ms"] = frame["stdev_s"] * 1.0e3
    return frame
