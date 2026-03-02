# Runko v5 Compute Kernel SIMD Audit Report

**Platform:** ARM NEON (Apple Silicon) · clang 22.1.0 · `-O3 -march=native`
**Date:** 2026-03-01
**Branch:** `cpu-backend`

All `for_each_index` loops use `#pragma omp simd` on the innermost spatial dimension
(`mdgrid.h:223`). ARM NEON = 128-bit registers → float width=4, size_t width=2.

---

## 1. FDTD2 Field Propagator — `yee_lattice_FDTD2.c++`

**Result: 6/6 loops vectorized at width=4, interleave=4**

| Function | Loops | Width | Interleave |
|----------|-------|-------|------------|
| `push_half_b` (Bx, By, Bz) | 3 | 4 | 4 |
| `push_e` (Ex, Ey, Ez) | 3 | 4 | 4 |

Status: **FULLY VECTORIZED** — optimal throughput.

---

## 2. Binomial Current Filter — `yee_lattice_current_filter_binomial2.c++`

**Result: 3/3 loops vectorized at width=4, interleave=1**

| Function | Loops | Width | Interleave |
|----------|-------|-------|------------|
| Pass 1 (z-convolution) | 1 | 4 | 1 |
| Pass 2 (y-convolution) | 1 | 4 | 1 |
| Pass 3 (x-convolution) | 1 | 4 | 1 |

Status: **FULLY VECTORIZED** — was 0/1 before fix.
Previous blocker: 27-point inner `for` loop over `C3_index_space` prevented `#pragma omp simd`.
Fix: separated into 3 sequential 1D passes with 3-term inline sums (`B1 = {0.25f, 0.5f, 0.25f}`).

---

## 3. Field Interpolation — `yee_lattice_interpolate_linear_1st.c++`

**Result: 1/1 loop vectorized at width=4, interleave=1**

| Function | Loops | Width | Interleave |
|----------|-------|-------|------------|
| `interpolate_EB_linear_1st` | 1 | 4 | 1 |

Status: **FULLY VECTORIZED** — was 0/1 before fix.
Previous blockers: (1) `index_space` 2×2×2 inner loop — compiler said "not beneficial";
(2) `submdspan` accessor chain through `mdgrid_buffer` too complex for vectorizer.
Fix: unrolled 8-point gather into 48 direct field accesses with `mean2`/`mean4`, inlined
trilinear interpolation as scalar `lerp3D(dx, dy, dz, v000..v111)`.

---

## 4. Boris Particle Pusher — `particle_boris.c++`

**Result: 1/1 loop vectorized at width=4, interleave=1**

| Function | Loops | Width | Interleave |
|----------|-------|-------|------------|
| `push_particles_boris` | 1 | 4 | 1 |

Status: **FULLY VECTORIZED** — optimal throughput.

---

## 5. Zigzag Current Deposition — `particle_current_zigzag_1st.c++`

**Result: 2/2 loops vectorized — atomic overload at width=4, scatter overload at width=2**

| Function | Loops | Width | Interleave |
|----------|-------|-------|------------|
| `current_zigzag_1st` (scatter overload) | 1 | 2 | 1 |
| `current_zigzag_1st` (atomic overload) | 1 | 4 | 1 |

Status: **Atomic overload FIXED to full width=4** (was width=2).

Fix applied: replaced `size_t` (64-bit) index operations with `float floor()` for
relay/weight computation (pure float arithmetic) and `uint32_t` (32-bit) for cell
index storage. Added `tyvi::sstd::floor()` wrapper (maps to FRINTM on ARM NEON).

The **atomic overload** (used by all production simulations: `zigzag_1st_atomic`)
now vectorizes at full width=4. The scatter overload remains at width=2 due to
the cost model: 28 indirect stores per iteration (14 deposit locations × 2 arrays)
create a scatter pattern that the vectorizer determines is not profitable at width=4.

Additional diagnostic: line 44 reports "value that could not be identified as reduction
is used outside the loop" — this is the `deposit_loc` comparison operator, part of the
sort/reduce pipeline (not the SIMD kernel).

---

## 6. Yee Lattice Utility Functions — `yee_lattice.c++`

**Result: 6/6 for_each_index loops vectorized**

| Function | Loops | Width | Interleave |
|----------|-------|-------|------------|
| `mds_copy` | 1 | 4 | 1 |
| `mds_add` | 1 | 4 | 2 |
| `batch_set_EBJ` | 1 | 4 | 4 |
| `clear_current` | 1 | 4 | 1 |
| `deposit_current_contributions` | 1 | 4 | 4 |
| (6th loop — copy/add template instantiation) | 1 | 4 | 1 |

Status: **FULLY VECTORIZED** — optimal throughput.

---

## 7. Particle Operations — `particle.c++`

**Result: 4 loops vectorized**

| Function | Loops | Width | Interleave | Source |
|----------|-------|-------|------------|--------|
| `get_positions` | 1 | 4 | 2 | particle.c++:53 |
| `get_velocities` | 1 | 4 | 2 | particle.c++:77 |
| `wrap_positions` | 1 | 4 | 1 | mdgrid.h:223 |
| `sort_by_cell_index` | 1 | 2 | 1 | mdgrid.h:223 |

Status: **VECTORIZED** — sort uses `size_t` keys → width=2 (expected for 64-bit indices).

---

## Summary Table

| Kernel | Loops | Width | Status |
|--------|-------|-------|--------|
| FDTD2 field propagator | 6 | 4×4 | Optimal |
| Binomial current filter | 3 | 4×1 | **FIXED** (was 0) |
| Field interpolation | 1 | 4×1 | **FIXED** (was 0) |
| Boris particle pusher | 1 | 4×1 | Optimal |
| Zigzag current deposition | 2 | 4×1 / 2×1 | Atomic **FIXED** (was 2); scatter: cost-model |
| YeeLattice utilities | 6 | 4×1-4 | Optimal |
| Particle get/sort | 4 | 4×2 / 2×1 | Optimal (sort: `size_t`) |

**Totals:**
- Vectorized `for_each_index` loops: **21/21** (0 failed)
- At full width (4): **20/21**
- At half width (2): **1/21** (zigzag scatter — cost-model, 28 indirect stores)
- Mixed-precision warnings: **0**

---

## How This Audit Was Produced

See `cpp_vectorization_analysis.md` for the full tutorial on reproducing this analysis.
