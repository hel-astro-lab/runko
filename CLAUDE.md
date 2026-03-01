# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Runko v5 is a modern C++23/Python plasma simulation framework for astrophysics. High-performance C++ kernels are exposed to Python via pybind11, enabling rapid prototyping with efficient execution. The grid infrastructure uses the [corgi](https://github.com/natj/corgi) library for tile-based domain decomposition with MPI parallelization (tested up to ~20k cores).

## Python Environment

A Python venv must be active before building or running runko. Activate it with:

```bash
source venv/runko5/bin/activate
```

The venv includes: numpy, scipy, matplotlib, mpi4py (built against the local OpenMPI / Homebrew clang).

## Build Commands

The CPU backend uses `mpic++` as the C++ compiler, which wraps Homebrew clang 22.1.0 (`/opt/homebrew/opt/llvm/bin/clang++`) with OpenMPI (`~/bin/openmpi51/`).

```bash
# Configure using CMake presets (preferred)
cmake --preset unixlike-clang-cpu-debug     # CPU debug build
cmake --preset unixlike-clang-cpu-release   # CPU release build

# Or configure manually (from repo root)
cmake -B build -DCMAKE_CXX_COMPILER=mpic++ -DCMAKE_BUILD_TYPE=Debug -DTYVI_BACKEND=cpu

# Build (compiles C++ and produces runko_cpp_bindings.so in repo root)
cmake --build build

# Run single-rank unit tests
python -m unittest discover -s tests/ -v

# Run a single test file
python -m unittest tests/test_emf.py -v

# Run multi-rank MPI tests (via CTest)
cd build && ctest -V

# Run a specific multi-rank test
mpirun -np 4 python tests/multirank_tests/test_pic_simulation.py
```

The shared library `runko_cpp_bindings*.so` is output to the repository root directory. Tests and simulations must be run from a location where Python can find it.

## Architecture

### Hybrid C++/Python Design

- **C++ kernels** (`core/`, `tools/`, `io/`): Performance-critical physics solvers, data structures, and I/O
- **pybind11 bindings** (`bindings/`): Single shared library `runko_cpp_bindings` exposing all C++ modules to Python
- **Python package** (`runko/`): High-level simulation API (Simulation, TileGrid, Configuration classes)
- **Projects** (`projects/`): User simulation scripts that drive the framework

### Core Physics Modules (`core/`)

- **`emf/`** — 3D FDTD electromagnetic field solver on a staggered Yee lattice. Includes propagator (FDTD2), current filter (binomial2), and field interpolation (linear 1st order).
- **`pic/`** — 3D3V Particle-In-Cell module. Boris particle pusher, zigzag current deposition, MPI particle communication.
- **`vlv/`**, **`ffe/`**, **`qed/`** — Vlasov, force-free electrodynamics, and QED radiation modules exist but are currently commented out in `bindings/CMakeLists.txt`.

### Grid Infrastructure

The simulation domain is decomposed into **tiles** managed by the `corgi` library. Each tile contains local field data (Yee lattice) and particle containers. Tiles have a halo region (size=3) for neighbor communication. Domain partitioning uses either Hilbert curve or caterpillar track ordering.

### Simulation Flow

A typical simulation script (see `projects/pic-turbulence/pic.py`):
1. Create `Configuration` with grid and physics parameters
2. Construct `TileGrid` — distributes tiles across MPI ranks
3. Initialize tiles with field/particle data via `batch_set_EBJ` and `batch_inject_to_cells`
4. Call `tile_grid.configure_simulation(config)` — sets up halos, virtual tiles, I/O
5. Define a lap function using the `MethodWrapper` API with prefixed commands:
   - `grid_*` — field operations (push_half_b, push_e, filter_current, etc.)
   - `prtcl_*` — particle operations (push, sort, deposit_current)
   - `comm_*` — MPI communication (local, external) with `comm_mode` flags
   - `io_*` — output (emf_snapshot, average_kinetic_energy, etc.)
6. Run with `simulation.for_each_lap(lap_function)`

### External Dependencies (git submodules in `external/`)

- **corgi** — Tile-based grid with MPI load-balancing
- **tyvi** — mdspan-based multi-dimensional grid buffers
- **cppitertools**, **ezh5**, **fmt** — Utility libraries

System dependencies: MPI, HDF5, Python (with mpi4py, numpy, pybind11), OpenMP.

## Code Conventions

- C++23 standard with modern features (concepts, ranges, mdspan)
- Formatting: LLVM-based `.clang-format` with 88-column limit, 2-space indent, WebKit brace style, no tabs
- C++ source files use `.c++` extension (not `.cpp`)
- Single-precision floats by default; mixed-precision for critical numerical parts
- Python tests use `unittest` framework
- All physics is 3D (the `threeD` namespace appears throughout bindings and Python imports)

### Naming Conventions

- **Variables and functions**: `snake_case` — e.g., `ex`, `gam`, `get_container()`
- **Classes and objects**: `CamelCase` with leading capital — e.g., `Tile`, `ParticleContainer`
- Prefer short but descriptive names that mirror the underlying mathematical formulas
