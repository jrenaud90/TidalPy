"""Run the full performance suite and append a timestamped JSON record under ``results/``.

Run it from this directory so the installed TidalPy is imported rather than the source tree:

    cd Benchmarks_x/Performance
    python run_all.py
    python run_all.py --label "post-eos-cache"

Each run captures the current environment (date, TidalPy version, git commit, OS, CPU, Python) so
the timings can later be tracked across all of those axes in ``Perf_Trends.ipynb``.
"""
from __future__ import annotations

import argparse

import harness
import tasks  # noqa: F401  (importing registers the tasks through the @benchmark decorator)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the TidalPy performance suite.")
    parser.add_argument("--label", default="", help="Optional label stored with this run.")
    args = parser.parse_args()

    task_list = harness.registered_tasks()
    if not task_list:
        print(
            "No benchmark tasks are registered yet. Add tasks in "
            "tasks/common_tasks.py, then rerun."
        )
        return 0

    print(f"Running {len(task_list)} performance task(s)...")
    path = harness.run_suite(label=args.label)
    print(f"Done. Results appended at {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
