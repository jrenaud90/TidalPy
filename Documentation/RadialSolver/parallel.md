# Parallel Love Solves

_Updated: 2026-09-28_

Both entry points to the radial solver, the standalone `radial_solver` and `LayeredWorld.solve_love_numbers`, release Python's interpreter lock while the solve runs, so a thread pool can run several solves at once. A process pool works as well, with the setup described below. `calc_tides` and the 3D grid methods also use threads of their own, described at the end.

The speedups quoted below were measured on a 16-thread desktop with 64 forcing frequencies and 16 threads. They depend on the problem size and the machine, so time your own workload.

## Standalone `radial_solver`

Each `radial_solver` call builds its own temporary world from the arrays it is given, so two calls share no state and can run on any number of threads:

```python
import os
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from TidalPy.RadialSolver import build_rs_input_homogeneous_layers, radial_solver
from TidalPy.Rheology import Elastic, Maxwell


def solve_k2(frequency):
    """Build the inputs for one forcing frequency and return the tidal k2."""
    # Use a helper function to build our radial solver inputs (see "build_inputs.md" for details on these helpers)
    build_data = build_rs_input_homogeneous_layers(
        6000.0e3,                                         # Planet radius [m]
        frequency,                                        # Forcing frequency [rad s-1]
        density_tuple=(8000.0, 3300.0),
        static_bulk_modulus_tuple=(2.0e11, 1.0e11),
        static_shear_modulus_tuple=(1.0e11, 5.0e10),
        bulk_viscosity_tuple=(1.0e18, 1.0e18),
        shear_viscosity_tuple=(1.0e20, 1.0e19),
        layer_type_tuple=("solid", "solid"),
        layer_is_static_tuple=(False, False),
        layer_is_incompressible_tuple=(False, False),
        shear_rheology_model_tuple=Maxwell(),
        bulk_rheology_model_tuple=Elastic(),
        radius_fraction_tuple=(0.5, 1.0),
        slices_tuple=(20, 30))

    solution = radial_solver(
        *build_data,
        degree_l=2,
        raise_on_fail=True)   # A failed solve raises in the calling thread
    return solution.k


frequencies = np.logspace(-6, -3, 64)   # [rad s-1]
with ThreadPoolExecutor(max_workers=os.cpu_count()) as pool:
    k2 = list(pool.map(solve_k2, frequencies))   # In the order of `frequencies`
```

The argument checks and the construction of the returned solution run in Python with the interpreter lock, so small solves gain less than large ones. This 50-slice problem ran 2.5 times faster on 16 threads, and a 300-slice problem using `integration_rtol=1e-10` ran ~5 times faster.

## World-Attached Solves

A world takes one of its heavy calls at a time. `solve_eos`, `solve_love_numbers`, `calc_tides`, and the 3D calls on one world wait for each other, while calls on separate worlds run in parallel (see the thread notes in [Worlds](../Structures/worlds/worlds.md)). Threads sharing one world therefore gain no speed, and they also need care to read the right result.

### Race Conditions

A world keeps the result of its most recent solve, and `solve_love_numbers` returns its dict by reading that stored result after the solve finishes. The world's lock covers the solve and each read on its own, not the pair. When several threads call `solve_love_numbers` on the same world, another thread's solve can replace the stored result between one thread's solve and its read, and that thread returns Love numbers for the wrong frequency with no error or warning. In a test of 64 frequencies on 16 threads sharing the bundled Io, 1219 of 1280 returned results belonged to another thread's solve.

The same applies to every solve followed by a read of what it stored:

- `solve_love_numbers` followed by `love_number_k`, `get_love_radial_y`, or the other Love-number getters.
- `calc_tides` followed by `get_tidal_heating`, `get_tidal_love_k`, or the other tidal getters.
- `solve_eos` followed by the profile getters (`get_density`, `get_gravity`, and the rest).

Changing a world's configuration (`set_tide_config`, `set_tide_model`, a layer's models, or `set_spin_frequency`) from one thread while another thread solves on it is also a race.

There are two ways to avoid these races: give each thread its own world, or hold a Python lock around each solve and the reads of its result.

### One World per Thread

Separate worlds share no state, so this is the approach that also runs in parallel. Build the worlds before starting the threads, and give each thread one world:

```python
import os
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from TidalPy.Structures import build_world

num_threads = os.cpu_count()
frequencies = np.logspace(-6, -3, 64)

worlds = []
for _ in range(num_threads):
    world = build_world("io")   # One independent world per thread
    world.solve_eos()           # Each world solves its own interior
    worlds.append(world)


def solve_chunk(thread_index):
    """Solve every frequency assigned to this thread on this thread's own world."""
    world = worlds[thread_index]
    return [world.solve_love_numbers(frequency=frequency)["love_number_k"]
            for frequency in frequencies[thread_index::num_threads]]


with ThreadPoolExecutor(max_workers=num_threads) as pool:
    chunks = list(pool.map(solve_chunk, range(num_threads)))

k2 = np.empty(frequencies.size, dtype=np.complex128)
for thread_index, chunk in enumerate(chunks):
    k2[thread_index::num_threads] = chunk                 # Back into the order of `frequencies`
```

This ran 6.6 times faster than a serial loop on one world, with identical results. To copy a world that was built or edited in code rather than loaded by name, rebuild it from its configuration with `build_world(world.get_config_dict())`. The copy starts unsolved, so you need to call `solve_eos` on it again.

### Sharing One World

A `threading.Lock` held around each solve and every read of its result also removes the race, with one world shared by every thread.

```python
import os
import threading
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from TidalPy.Structures import build_world

frequencies = np.logspace(-6, -3, 64)   # [rad s-1]

io = build_world("io")
io.solve_eos()
io_lock = threading.Lock()


def solve_k2_shared(frequency):
    """Solve and read under one lock, so no other thread's solve can come between them."""
    with io_lock:
        result = io.solve_love_numbers(frequency=frequency)
        surface_y1 = io.get_love_radial_y(
            io.radius,
            y_idx=0)   # Other reads of this solve go inside the lock too
    return result["love_number_k"], surface_y1


with ThreadPoolExecutor(max_workers=os.cpu_count()) as pool:
    results = list(pool.map(solve_k2_shared, frequencies))   # In the order of `frequencies`
```

The lock lets only one solve run at a time, so the threads wait on each other instead of working in parallel. These 64 solves took 25.9 ms on 16 threads against 22.5 ms for a serial loop, 15 percent slower. One world per thread took 3.4 ms. Building and solving another bundled Io took about 0.8 ms, so extra worlds are usually the better choice. The lock suits a world too expensive to copy, or code where other threads must read one world's results.

## Process Pools

A process pool also runs world-attached and standalone solves in parallel, and each process has its own worlds and configuration, so none of the races above can occur between processes. Three points differ from threads:

- Worlds cannot be pickled, so they cannot be sent to a worker. Send what the world is built from instead: a bundled name, a TOML path, or the dict from `world.get_config_dict()`, and build and solve the world inside the worker.
- On Windows and macOS, workers are started fresh (the `spawn` method) and load the configuration file, so a `TidalPy.reinit(provided_config=...)` override made in the main process does not reach them. Pass the same override to each worker through the pool's `initializer`.
- The script's main code must sit under `if __name__ == "__main__":`, because each spawned worker imports the script.

```python
from concurrent.futures import ProcessPoolExecutor

import numpy as np

import TidalPy
from TidalPy.Structures import build_world


def init_worker(config_overrides):
    """Apply the main process's configuration overrides once in each worker."""
    TidalPy.reinit(provided_config=config_overrides)


def solve_chunk(world_source, frequency_chunk):
    """Build and solve a world inside the worker, then solve each assigned frequency."""
    world = build_world(world_source)
    world.solve_eos()
    return [world.solve_love_numbers(frequency=frequency)["love_number_k"] for frequency in frequency_chunk]


if __name__ == "__main__":
    config_overrides = {
        "radial_solver": {"rtol": 1.0e-8},
        "numerical": {"love_solve_threads": 1}}           # The pool's workers already occupy the machine
    TidalPy.reinit(provided_config=config_overrides)      # Applies to this process only
    world_source = build_world("io").get_config_dict()    # A bundled name or a TOML path also works

    num_workers = 8
    frequencies = np.logspace(-6, -3, 64)                 # [rad s-1]
    chunks = [frequencies[worker::num_workers] for worker in range(num_workers)]
    with ProcessPoolExecutor(
            max_workers=num_workers,
            initializer=init_worker,
            initargs=(config_overrides,)) as pool:
        results = list(pool.map(solve_chunk, [world_source] * num_workers, chunks))
```

Starting a worker, importing TidalPy, and building its world take a few hundred milliseconds, which is more than the 64 solves above take on one thread. Processes can work for long jobs, such as a parameter sweep of thousands of solves or of full thermal-orbital evolutions. Threads are better for short ones.

## Threads Inside `calc_tides`

`calc_tides` runs one Love solve per unique (degree, frequency) pair of the tidal modes. With the radial-solver Love methods, once there are at least `love_solve_min_parallel` solves (default 3), it spreads them over `love_solve_threads` threads, both in the `[numerical]` section of the configuration (see [TidalPy Configurations](../Overview/2_TidalPy_Configurations.md)). The default, 0, uses `max(num_logical_processors - 4, 1)`. `love_solve_threads = 1` keeps every solve on a single thread. The results are identical for any thread count. The quasi-homogeneous Love methods take microseconds per solve and always runs on one thread.

The table times `calc_tides` on the bundled Io at an eccentricity of 0.05 and an obliquity of 0.1 [rad], with 1 and 16 Love-solve threads. The per-layer heating (`layer_tidal_heating`, on by default) then integrates the radial solutions over each layer on the same threads, at least `tides_3d_min_radii_per_thread` radii (default 8) to a thread. Part of that integral runs on the calling thread, so it gains less:

| Tide settings | Love solves | 1 thread | Change, 16 threads | Change without the per-layer heating |
|---|---|---|---|---|
| Degree 2, e^10, synchronous | 5 | 2.1 ms | 1.6x faster | 2.2x faster |
| Degree 2, e^10, obliquity level 2, spin 1.2 n | 28 | 12.7 ms | 3.1x faster | 5.5x faster |
| Degrees 2 to 3, e^10, obliquity level 2, spin 1.2 n | 58 | 28.1 ms | 3.8x faster | 7.1x faster |
| Degrees 2 to 4, e^20, obliquity level 4, spin 1.2 n | 181 | 97.9 ms | 4.3x faster | 8.4x faster |

The default value of `love_solve_min_parallel = 3` was born out of testing where two solves roughly broke even on a one-layer body (0.17 ms per solve). But 3+ had every case faster, which sets the default of `love_solve_min_parallel`.

The 3D grid methods (`calc_3d_tides`, `calc_3d_stress_strain`, `calc_3d_displacements`, and `get_3d_tidal_heating_array`) spread their radial solves (from `love_solve_min_parallel` of them on) and their per-point work the same way through their `num_threads` argument, whose default of 0 means the same count (see [3D Tidal Stress, Strain, and Heating](../Tides/multilayer_3d_heating.md)).

Inside a thread or process pool that already occupies the machine, these threads compete with the pool's workers. Set `love_solve_threads` to 1 in each worker (through the pool's `initializer`, as in the example above) and pass `num_threads=1` to the 3D methods. In a test with 16 worker processes each running `calc_tides`, leaving the automatic threads on cost about 3 percent.

## Configuration and Logging

The TidalPy configuration is shared by every thread in a process. Set the configuration before starting them, do not call `TidalPy.reinit` while other threads are solving.

The logger is safe to use from any thread (see [Logging](../Utilities/logging.md)). Warnings from solves on different threads can arrive in any order. `set_console_level("error")` from `TidalPy.Utilities.logging` hides the solver's conditioning warnings during a large parallel run.
