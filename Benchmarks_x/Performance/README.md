# Performance tracking

A small harness (standard library plus `psutil`) that times common TidalPy tasks
and records every run so performance can be followed as a function of **date**, **OS / machine**, and
**TidalPy version**.

## Files

```
Performance/
  harness.py         timing engine, environment capture, and result loading
  run_all.py         runs every registered task and writes one JSON record
  tasks/
    common_tasks.py  the tasks (mirror the operations taught in Demos_x/)
  results/           one JSON file per run (accumulates locally per machine; commit runs you want to keep)
  Perf_Trends.ipynb  loads results/*.json and plots timings vs date / version / OS
```

## Running a suite

Run it from this directory (not the repository root, so the installed TidalPy is imported rather
than the source tree):

```
cd Benchmarks_x/Performance
python run_all.py
python run_all.py --label "post-eos-cache"
```

Each run writes `results/<timestamp>_<host>_<version>.json`. The file holds one `environment` block
(timestamp, TidalPy version, git commit, OS, machine, CPU, Python) and one timing record per task
(auto-calibrated inner loop count, then per-call min / median / mean / standard deviation in seconds).

## Adding a task

Register a zero-argument callable in `tasks/common_tasks.py`:

```python
from harness import benchmark

@benchmark("build_layered_world", group="structures")
def _build_layered_world():
    LayeredWorld("planet", 6.371e6, 5.972e24)
```

Do one-time construction in a `setup=` callback or at module scope so the timed call measures the
operation itself, not its setup. Prefer tasks that mirror what the demos teach, so the tracked
numbers reflect real usage.

## Viewing trends

Open `Perf_Trends.ipynb`. It calls `harness.load_results()` to flatten every run into one row per
(run, task) and plots each task's per-call time against commit date, TidalPy version, and OS, so a
regression or a speedup shows up directly.
