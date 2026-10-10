# Parallel Love Solves

_Updated: 2026-10-02_

`radial_solver` and `BaseWorld.solve_love_numbers` release Python's interpreter lock while they solve, so a thread pool runs several solves at once. Process pools work too, with the setup below. `calc_tides` and the 3D grid methods use threads of their own. Speedups below were measured once on a 16-thread desktop, with 64 forcing frequencies and 16 threads.

## Standalone `radial_solver`

Each `radial_solver` call builds its own temporary world, so calls share no state and can run on any number of threads:

```python
import os
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from TidalPy.RadialSolver import build_rs_input_homogeneous_layers, radial_solver
from TidalPy.Rheology import Elastic, Maxwell


def solve_k2(frequency):
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

Argument checks and building the returned solution hold the interpreter lock, so small solves gain less: this 50-slice problem ran 2.5 times faster on 16 threads, a 300-slice one at `integration_rtol=1e-10` about 5 times.

## World-Attached Solves

A world runs one heavy call at a time (`solve_eos`, `solve_love_numbers`, `calc_tides`, and the 3D calls wait for each other), while separate worlds run in parallel (see [Worlds](../Structures/worlds/worlds.md)). Threads sharing one world gain no speed and can read the wrong result.

### Race Conditions

A world keeps its latest result, and `solve_love_numbers` reads that stored result after solving. If another thread's solve replaces it in between, a thread returns Love numbers for the wrong frequency with no error or warning: with 64 frequencies on 16 threads sharing Io, 1219 of 1280 results belonged to another thread's solve. The same holds for every solve followed by a read of what it stored:

- `solve_love_numbers`, then the Love-number getters (`love_number_k`, `get_love_radial_y`, ...).
- `calc_tides`, then the tidal getters (`get_tidal_heating`, `get_tidal_love_k`, ...).
- `solve_eos`, then the profile getters (`get_density`, `get_gravity`, ...).

Changing a world's configuration (`set_tide_config`, `set_tide_model`, a layer's models, `set_spin_frequency`) while another thread solves on it is also a race. Give each thread its own world, or hold a lock around each solve and its reads.

### One World per Thread

Build the worlds before starting the threads and give each thread one:

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
    world.solve_eos()
    worlds.append(world)


def solve_chunk(thread_index):
    world = worlds[thread_index]
    return [world.solve_love_numbers(frequency=frequency)["love_number_k"]
            for frequency in frequencies[thread_index::num_threads]]


with ThreadPoolExecutor(max_workers=num_threads) as pool:
    chunks = list(pool.map(solve_chunk, range(num_threads)))

k2 = np.empty(frequencies.size, dtype=np.complex128)
for thread_index, chunk in enumerate(chunks):
    k2[thread_index::num_threads] = chunk     # Back into the order of `frequencies`
```

This ran 6.6 times faster than a serial loop on one world, with identical results. To copy a world built or edited in code, use `world.copy()` (or `copy.deepcopy(world)`), which carries every layer, model, and setting. The copy starts unsolved, so call `solve_eos` on it.

### Sharing One World

A `threading.Lock` around each solve and its reads also removes the race, with one world shared by every thread:

```python
import threading

from TidalPy.Structures import build_world

io = build_world("io")
io.solve_eos()
io_lock = threading.Lock()


def solve_k2_shared(frequency):
    with io_lock:   # The solve and every read of its result go inside the lock
        result = io.solve_love_numbers(frequency=frequency)
        surface_y1 = io.get_love_radial_y(io.radius, y_idx=0)
    return result["love_number_k"], surface_y1
```

The lock runs one solve at a time: these 64 solves took 25.9 ms on 16 threads, 22.5 ms serially, and 3.4 ms with one world per thread. Another Io takes about 0.8 ms to build and solve, so extra worlds are usually better; a lock suits worlds that are very expensive to build or copy.

## Process Pools

Each process has its own worlds and configuration, so the races above cannot occur between processes. Three points differ from threads:

- A world pickles and arrives in the worker unsolved: call `solve_eos` there. Sending what the world is built from (a bundled name, a TOML path, or `world.get_config_dict()`) works as well.
- On Windows and macOS, workers start fresh (the `spawn` method) and load the configuration file, so a `TidalPy.reinit(provided_config=...)` override in the main process does not reach them. Pass it to each worker through the pool's `initializer`.
- The script's main code must sit under `if __name__ == "__main__":`, because each spawned worker imports the script.

```python
from concurrent.futures import ProcessPoolExecutor

import numpy as np

import TidalPy
from TidalPy.Structures import build_world


def init_worker(config_overrides):
    TidalPy.reinit(provided_config=config_overrides)


def solve_chunk(world_source, frequency_chunk):
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

Starting a worker, importing TidalPy, and building its world take a few hundred milliseconds, more than the 64 solves above take on one thread. Use processes for long jobs (sweeps of thousands of solves, full thermal-orbital evolutions) and threads for short ones.

## Threads Inside `calc_tides`

`calc_tides` runs one Love solve per unique (degree, frequency) pair of the tidal modes. With the radial-solver Love methods and at least `love_solve_min_parallel` solves (default 3, where threading first came out faster), it spreads them over `love_solve_threads` threads (both `[numerical]` settings; see [TidalPy Configurations](../Overview/2_TidalPy_Configurations.md)). The default, 0, uses `max(num_logical_processors - 4, 1)`; 1 keeps every solve on one thread. Results are identical for any thread count. The quasi-homogeneous Love methods always run on one thread. The per-layer heating (`layer_tidal_heating`, on by default) uses the same threads, at least `tides_3d_min_radii_per_thread` radii (default 8) each, but runs partly on the calling thread, so it gains less.

On the bundled Io (eccentricity 0.05, obliquity 0.1 rad), 16 threads made `calc_tides` 1.6 times faster with 5 Love solves (2.1 ms on one thread) and 4.3 times faster with 181 (97.9 ms); without the per-layer heating, 2.2 and 8.4 times.

The 3D grid methods (`calc_3d_tides`, `calc_3d_stress_strain`, `calc_3d_displacements`, and `get_3d_tidal_heating_array`) spread their radial solves (from `love_solve_min_parallel` of them on) and their per-point work the same way through their `num_threads` argument, whose default of 0 means the same count (see [3D Tidal Stress, Strain, and Heating](../Tides/multilayer_3d_heating.md)).

Inside a thread or process pool that already occupies the machine, these threads compete with the pool's workers. Set `love_solve_threads` to 1 in each worker (through the pool's `initializer`, as above) and pass `num_threads=1` to the 3D methods. With 16 worker processes each running `calc_tides`, leaving the automatic threads on cost about 3 percent.

## Configuration and Logging

Every thread in a process shares the TidalPy configuration. Set it before starting the threads, and do not call `TidalPy.reinit` while other threads are solving.

The [logger](../Utilities/logging.md) is thread-safe, but warnings from different threads arrive in any order. `set_console_level("error")` from `TidalPy.Utilities.logging` hides the conditioning warnings during a large run.
