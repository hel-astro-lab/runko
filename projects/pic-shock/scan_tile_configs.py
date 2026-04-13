#!/usr/bin/env python3
"""
Scan optimal tile/mesh configurations for 3D PIC shock simulations.

Enumerates valid (Nx, Ny, Nz, NxMesh, NyMesh, NzMesh) factorizations,
computes memory, load-balance, and GPU efficiency metrics,
and ranks configurations by a weighted combined score.

Usage:
    python scan_tile_configs.py --Lx 65536 --Ly 128 --Lz 128 --n_gpus 8 --ppc 2
    python scan_tile_configs.py --Lx 65536 --Ly 128 --Lz 128 --n_gpus 8 --plot
"""

import argparse
import math
import numpy as np
from itertools import product

try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, **kwargs):
        return iterable

from runko.hilbert import Hilbert3D


# ============================================================================
# Section 2: Configuration enumeration
# ============================================================================

def enumerate_factorizations(L, min_mesh=8):
    """Return list of (N, NMesh) where N is a power of 2 and NMesh = L/N >= min_mesh."""
    results = []
    p = 1
    while p <= L:
        if L % p == 0:
            nmesh = L // p
            if nmesh >= min_mesh:
                results.append((p, nmesh))
        p *= 2
    return results


# ============================================================================
# Section 3: Memory model
# ============================================================================

HALO = 3
BYTES_PER_FIELD_CELL = 36   # 3 fields × 3 components × 4 bytes
BYTES_PER_PARTICLE   = 24   # (3 pos + 3 vel) × 4 bytes

def emf_bytes(nxm, nym, nzm):
    return BYTES_PER_FIELD_CELL * (nxm + 2*HALO) * (nym + 2*HALO) * (nzm + 2*HALO)

def virtual_emf_bytes(nxm, nym, nzm):
    """Memory for one hollow virtual tile: E + B (h=3) + J (h=6)."""
    h_eb = HALO
    shell_eb = nxm*nym*nzm - (nxm - 2*h_eb)*(nym - 2*h_eb)*(nzm - 2*h_eb)
    eb_bytes = 2 * 3 * shell_eb * 4  # E + B, 3 components, 4 bytes/float

    h_j  = 2 * HALO
    jx, jy, jz = nxm + 2*HALO, nym + 2*HALO, nzm + 2*HALO
    shell_j = jx*jy*jz - (jx - 2*h_j)*(jy - 2*h_j)*(jz - 2*h_j)
    j_bytes = 3 * shell_j * 4

    return eb_bytes + j_bytes

def particle_bytes(nxm, nym, nzm, ppc_eff, n_species):
    return n_species * ppc_eff * nxm * nym * nzm * BYTES_PER_PARTICLE


# ============================================================================
# Section 4: Partition and virtual tile counting
# ============================================================================

def build_partition(Nx, Ny, Nz, n_gpus):
    """Build Hilbert-curve partition grid; returns igrid[Nx, Ny, Nz] of rank owners."""
    m0 = int(math.log2(Nx))
    m1 = int(math.log2(Ny))
    m2 = int(math.log2(Nz))

    total_tiles = Nx * Ny * Nz
    hgen = Hilbert3D(m0, m1, m2)

    idx = np.arange(total_tiles, dtype=np.int64)
    x = idx // (Ny * Nz)
    yz = idx % (Ny * Nz)
    y = yz // Nz
    z = yz % Nz

    grid = hgen.hindex_vec(x, y, z).reshape(Nx, Ny, Nz)

    # hmax = total_tiles - 1 (verified analytically for general Hilbert curves)
    igrid = (n_gpus * grid) // total_tiles
    return igrid.astype(int)


MOORE_OFFSETS = [
    (di, dj, dk)
    for di in (-1, 0, 1)
    for dj in (-1, 0, 1)
    for dk in (-1, 0, 1)
    if not (di == 0 and dj == 0 and dk == 0)
]

def count_virtual_tiles(igrid, Nx, Ny, Nz, n_gpus):
    """Return array of virtual tile counts per rank (periodic wrapping on all axes)."""
    vtiles = np.zeros(n_gpus, dtype=int)

    for rank in range(n_gpus):
        owned = (igrid == rank)
        dilated = owned.copy()
        for di, dj, dk in MOORE_OFFSETS:
            shifted = np.roll(owned, di, axis=0)
            shifted = np.roll(shifted, dj, axis=1)
            shifted = np.roll(shifted, dk, axis=2)
            dilated |= shifted
        vtiles[rank] = np.sum(dilated & ~owned)

    return vtiles


# ============================================================================
# Section 5: Load balance scoring
# ============================================================================

def load_balance_score(igrid, Nx, n_gpus, va_over_c):
    """Score load balance across shock evolution; returns (score, worst_imbalance)."""
    worst_imbalance = 1.0

    for f_inj in np.linspace(0.05, 1.0, 20):
        ix_inj   = min(int(np.ceil(f_inj * Nx)), Nx)
        ix_shock = min(int(np.ceil(f_inj * va_over_c * Nx)), ix_inj)

        loads = np.zeros(n_gpus)
        for rank in range(n_gpus):
            # downstream tiles (ix < ix_shock): weight 4
            if ix_shock > 0:
                down_mask = igrid[:ix_shock, :, :] == rank
                loads[rank] += 4.0 * np.sum(down_mask)
            # upstream tiles (ix_shock <= ix < ix_inj): weight 1
            if ix_inj > ix_shock:
                up_mask = igrid[ix_shock:ix_inj, :, :] == rank
                loads[rank] += 1.0 * np.sum(up_mask)

        total_load = np.sum(loads)
        if total_load == 0:
            continue

        ideal = total_load / n_gpus
        imbalance = np.max(loads) / ideal
        worst_imbalance = max(worst_imbalance, imbalance)

    score = 1.0 / worst_imbalance
    return score, worst_imbalance


# ============================================================================
# Section 6: GPU efficiency scoring
# ============================================================================

def halo_efficiency(nxm, nym, nzm):
    """Fraction of allocated cells that are interior (non-halo). 1.0 = no waste."""
    interior = nxm * nym * nzm
    total    = (nxm + 2*HALO) * (nym + 2*HALO) * (nzm + 2*HALO)
    return interior / total


def aspect_ratio_score(nxm, nym, nzm):
    """Penalize extremely elongated tiles. 1.0 = cubic, 0.0 = degenerate."""
    dims = sorted([nxm, nym, nzm])
    # ratio of geometric mean to max dimension
    gmean = (dims[0] * dims[1] * dims[2]) ** (1.0/3.0)
    return gmean / dims[2]


def gpu_efficiency_score(nxm, nym, nzm, tiles_per_gpu):
    cells_per_tile = nxm * nym * nzm
    size_score   = min(1.0, cells_per_tile / 32768.0)
    launch_score = min(1.0, 4.0 / max(tiles_per_gpu, 1)) if tiles_per_gpu > 4 else 1.0
    halo_eff     = halo_efficiency(nxm, nym, nzm)
    return size_score * launch_score * halo_eff


# ============================================================================
# Section 7: Combined scoring and output
# ============================================================================

def scan_configs(args):
    Lx, Ly, Lz = args.Lx, args.Ly, args.Lz
    n_gpus      = args.n_gpus
    ppc         = args.ppc
    n_species   = args.n_species
    sigma       = args.sigma
    gpu_mem_lim = args.gpu_mem_gb * 1e9
    min_mesh    = args.min_mesh
    w_mem, w_bal, w_gpu = args.w_mem, args.w_bal, args.w_gpu

    va_over_c = math.sqrt(sigma / (1.0 + sigma))
    f_down = va_over_c  # fraction of filled domain that is downstream
    print(f"sigma = {sigma:.1f}  =>  v_A/c = {va_over_c:.4f}  =>  f_downstream = {f_down:.3f}")
    print(f"Domain: {Lx} x {Ly} x {Lz} cells,  n_gpus = {n_gpus},  ppc = {ppc}")
    print()

    fx = enumerate_factorizations(Lx, min_mesh)
    fy = enumerate_factorizations(Ly, min_mesh)
    fz = enumerate_factorizations(Lz, min_mesh)

    print(f"Factorizations: x:{len(fx)}  y:{len(fy)}  z:{len(fz)}  total: {len(fx)*len(fy)*len(fz)}")

    results = []
    n_skipped_tiles = 0
    n_skipped_mem   = 0

    # y and z are symmetric; keep only (Ny,NyM) <= (Nz,NzM) to avoid duplicates
    configs = [(a, b, c) for a, b, c in product(fx, fy, fz) if b <= c]
    for (nx, nxm), (ny, nym), (nz, nzm) in tqdm(configs, desc="Scanning configs"):
        total_tiles = nx * ny * nz
        tiles_per_gpu = total_tiles / n_gpus

        # need at least 1 tile per GPU
        if total_tiles < n_gpus:
            n_skipped_tiles += 1
            continue

        # memory per tile
        emf_per_tile  = emf_bytes(nxm, nym, nzm)
        vemf_per_tile = virtual_emf_bytes(nxm, nym, nzm)

        # worst-case particle memory: weighted avg of downstream + upstream
        # at steady state: f_down fraction downstream (4x), rest upstream (1x)
        ppc_down = 4.0 * ppc
        ppc_up   = 1.0 * ppc
        prtcl_down = particle_bytes(nxm, nym, nzm, ppc_down, n_species)
        prtcl_up   = particle_bytes(nxm, nym, nzm, ppc_up,   n_species)

        # early skip: lower-bound memory assuming perfect balance, zero virtual tiles
        mem_lower_bound = tiles_per_gpu * (
            emf_per_tile + f_down * prtcl_down + (1 - f_down) * prtcl_up
        )
        if mem_lower_bound > gpu_mem_lim:
            n_skipped_mem += 1
            continue

        # build Hilbert partition and count virtual tiles
        igrid = build_partition(nx, ny, nz, n_gpus)
        vtiles = count_virtual_tiles(igrid, nx, ny, nz, n_gpus)

        # per-rank memory: real tiles (EMF + particles) + virtual tiles (EMF only)
        tiles_per_rank = np.bincount(igrid.flatten(), minlength=n_gpus)

        # for each rank, count downstream vs upstream tiles (steady-state: full domain)
        ix_shock_ss = min(int(np.ceil(f_down * nx)), nx)
        rank_mem = np.zeros(n_gpus)
        for rank in range(n_gpus):
            n_down = np.sum(igrid[:ix_shock_ss, :, :] == rank)
            n_up   = tiles_per_rank[rank] - n_down
            rank_mem[rank] = (
                tiles_per_rank[rank] * emf_per_tile
                + n_down * prtcl_down
                + n_up   * prtcl_up
                + vtiles[rank] * vemf_per_tile
            )

        max_gpu_mem = np.max(rank_mem)

        if max_gpu_mem > gpu_mem_lim:
            n_skipped_mem += 1
            continue

        # scores
        mem_score = 1.0 - max_gpu_mem / gpu_mem_lim
        bal_score, worst_imbal = load_balance_score(igrid, nx, n_gpus, va_over_c)
        eff_score = gpu_efficiency_score(nxm, nym, nzm, tiles_per_gpu)

        halo_eff = halo_efficiency(nxm, nym, nzm)
        aspect   = aspect_ratio_score(nxm, nym, nzm)

        # weighted geometric mean: any near-zero component kills the score
        eps = 1e-6
        combined = (max(mem_score, eps)**w_mem
                  * max(bal_score, eps)**w_bal
                  * max(eff_score, eps)**w_gpu)

        results.append({
            'Nx': nx, 'Ny': ny, 'Nz': nz,
            'NxMesh': nxm, 'NyMesh': nym, 'NzMesh': nzm,
            'total_tiles': total_tiles,
            'tiles_per_gpu': tiles_per_gpu,
            'max_vtiles': int(np.max(vtiles)),
            'mean_vtiles': float(np.mean(vtiles)),
            'emf_MB': emf_per_tile / 1e6,
            'vemf_MB': vemf_per_tile / 1e6,
            'prtcl_down_MB': prtcl_down / 1e6,
            'prtcl_up_MB': prtcl_up / 1e6,
            'max_gpu_GB': max_gpu_mem / 1e9,
            'mem_score': mem_score,
            'bal_score': bal_score,
            'worst_imbal': worst_imbal,
            'halo_eff': halo_eff,
            'aspect': aspect,
            'eff_score': eff_score,
            'combined': combined,
        })

    print(f"Evaluated: {len(results)},  skipped (too few tiles): {n_skipped_tiles},  skipped (OOM): {n_skipped_mem}")
    print()

    # sort by combined score
    results.sort(key=lambda r: r['combined'], reverse=True)
    return results


# ============================================================================
# Section 8: Output formatting
# ============================================================================

def print_table(results, top=20):
    top_results = results[:top]

    hdr = (
        f"{'#':>3s} | {'Nx':>5s} {'Ny':>4s} {'Nz':>4s} | "
        f"{'NxM':>5s} {'NyM':>5s} {'NzM':>5s} | "
        f"{'T/GPU':>6s} {'Vt/GPU':>7s} | "
        f"{'EMF':>7s} {'Prt':>7s} | "
        f"{'GPU':>7s} | "
        f"{'Bal':>5s} {'Halo':>5s} {'Eff':>5s} | "
        f"{'Score':>6s}"
    )
    sep = "-" * len(hdr)
    units = (
        f"{'':>3s} | {'':>5s} {'':>4s} {'':>4s} | "
        f"{'':>5s} {'':>5s} {'':>5s} | "
        f"{'':>6s} {'':>7s} | "
        f"{'MB':>7s} {'MB':>7s} | "
        f"{'GB':>7s} | "
        f"{'':>5s} {'':>5s} {'':>5s} | "
        f"{'':>6s}"
    )

    print(sep)
    print(hdr)
    print(units)
    print(sep)

    for i, r in enumerate(top_results, 1):
        print(
            f"{i:3d} | "
            f"{r['Nx']:5d} {r['Ny']:4d} {r['Nz']:4d} | "
            f"{r['NxMesh']:5d} {r['NyMesh']:5d} {r['NzMesh']:5d} | "
            f"{r['tiles_per_gpu']:6.0f} {r['max_vtiles']:7d} | "
            f"{r['emf_MB']:7.2f} {r['prtcl_down_MB']:7.2f} | "
            f"{r['max_gpu_GB']:7.3f} | "
            f"{r['bal_score']:5.3f} {r['halo_eff']:5.3f} {r['eff_score']:5.3f} | "
            f"{r['combined']:6.4f}"
        )

    print(sep)


# ============================================================================
# Section 9: Optional plotting
# ============================================================================

def plot_results(results, gpu_mem_lim_gb):
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipping plots")
        return

    if len(results) < 2:
        print("Too few configurations to plot")
        return

    tpg   = np.array([r['tiles_per_gpu'] for r in results])
    gmem  = np.array([r['max_gpu_GB']    for r in results])
    score = np.array([r['combined']       for r in results])
    bal   = np.array([r['bal_score']      for r in results])
    eff   = np.array([r['eff_score']      for r in results])
    halo  = np.array([r['halo_eff']       for r in results])
    vtpg  = np.array([r['max_vtiles']     for r in results])
    mem_s = np.array([r['mem_score']      for r in results])

    # shared x-axis range (tiles/GPU, log scale)
    xmin, xmax = 0.5, 2.0 * np.max(tpg)

    scat_kw = dict(s=30, edgecolors='k', linewidths=0.3)

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))

    # ---- row 1 ----

    # panel (0,0): tiles/GPU vs GPU memory, colored by combined score
    sc = axes[0, 0].scatter(tpg, gmem, c=score, cmap='viridis', **scat_kw)
    axes[0, 0].set_xlabel('tiles / GPU')
    axes[0, 0].set_ylabel('max GPU memory (GB)')
    axes[0, 0].set_title('memory vs tile count')
    axes[0, 0].set_xscale('log')
    axes[0, 0].set_xlim(xmin, xmax)
    axes[0, 0].set_ylim(0, gpu_mem_lim_gb)
    plt.colorbar(sc, ax=axes[0, 0], label='combined score')

    # panel (0,1): tiles/GPU vs balance score, colored by GPU efficiency
    sc2 = axes[0, 1].scatter(tpg, bal, c=eff, cmap='plasma', **scat_kw)
    axes[0, 1].set_xlabel('tiles / GPU')
    axes[0, 1].set_ylabel('balance score')
    axes[0, 1].set_title('load balance vs tile count')
    axes[0, 1].set_xscale('log')
    axes[0, 1].set_xlim(xmin, xmax)
    axes[0, 1].set_ylim(-0.05, 1.05)
    plt.colorbar(sc2, ax=axes[0, 1], label='GPU efficiency')

    # panel (0,2): tiles/GPU vs combined score
    axes[0, 2].scatter(tpg, score, c='steelblue', **scat_kw)
    axes[0, 2].set_xlabel('tiles / GPU')
    axes[0, 2].set_ylabel('combined score')
    axes[0, 2].set_title('combined score vs tile count')
    axes[0, 2].set_xscale('log')
    axes[0, 2].set_xlim(xmin, xmax)
    axes[0, 2].set_ylim(-0.05, 1.05)

    # ---- row 2 ----

    # panel (1,0): tiles/GPU vs halo efficiency, colored by combined score
    sc3 = axes[1, 0].scatter(tpg, halo, c=score, cmap='viridis', **scat_kw)
    axes[1, 0].set_xlabel('tiles / GPU')
    axes[1, 0].set_ylabel('halo efficiency')
    axes[1, 0].set_title('halo overhead vs tile count')
    axes[1, 0].set_xscale('log')
    axes[1, 0].set_xlim(xmin, xmax)
    axes[1, 0].set_ylim(-0.05, 1.05)
    plt.colorbar(sc3, ax=axes[1, 0], label='combined score')

    # panel (1,1): tiles/GPU vs virtual tiles/GPU, colored by balance score
    sc4 = axes[1, 1].scatter(tpg, vtpg, c=bal, cmap='plasma', **scat_kw)
    axes[1, 1].set_xlabel('tiles / GPU')
    axes[1, 1].set_ylabel('max virtual tiles / GPU')
    axes[1, 1].set_title('communication overhead')
    axes[1, 1].set_xscale('log')
    axes[1, 1].set_xlim(xmin, xmax)
    axes[1, 1].set_yscale('log')
    axes[1, 1].set_ylim(0.5, 2.0 * np.max(vtpg))
    plt.colorbar(sc4, ax=axes[1, 1], label='balance score')

    # panel (1,2): balance score vs GPU efficiency, colored by combined score
    sc5 = axes[1, 2].scatter(bal, eff, c=score, cmap='viridis', **scat_kw)
    axes[1, 2].set_xlabel('balance score')
    axes[1, 2].set_ylabel('GPU efficiency')
    axes[1, 2].set_title('balance vs efficiency tradeoff')
    axes[1, 2].set_xlim(-0.05, 1.05)
    axes[1, 2].set_ylim(-0.05, 1.05)
    plt.colorbar(sc5, ax=axes[1, 2], label='combined score')

    # mark top-3 on key panels
    for ax in [axes[0, 0], axes[0, 2]]:
        for i in range(min(3, len(results))):
            r = results[i]
            lbl = f"{r['Nx']}x{r['Ny']}x{r['Nz']}"
            yval = r['max_gpu_GB'] if ax is axes[0, 0] else r['combined']
            ax.annotate(lbl, (r['tiles_per_gpu'], yval),
                        fontsize=7, ha='left', va='bottom')

    plt.tight_layout()
    fname = 'tile_config_scan.png'
    plt.savefig(fname, dpi=150)
    print(f"Plot saved to {fname}")
    plt.show()


# ============================================================================
# Section 9b: Robustness analysis via weight perturbation
# ============================================================================

def robustness_analysis(results, w_mem, w_bal, w_gpu, top, weight_var):
    """Re-score top configs over a Cartesian grid of perturbed weights.

    For each (w_mem, w_bal, w_gpu) sample the combined score, then report
    mean and std of the score distribution.  A small std means the ranking
    is robust to the exact choice of weights.
    """
    top_results = results[:top]
    if not top_results:
        return

    # build perturbation grid: 5 points per axis → 125 raw combos,
    # then normalize so weights sum to 1 and each lies in [0, 1]
    n_pts = 5
    lo, hi = 1.0 - weight_var, 1.0 + weight_var
    factors = np.linspace(lo, hi, n_pts)
    weight_combos = []
    for fm in factors:
        for fb in factors:
            for fg in factors:
                wm = max(0.0, w_mem * fm)
                wb = max(0.0, w_bal * fb)
                wg = max(0.0, w_gpu * fg)
                wsum = wm + wb + wg
                if wsum > 0:
                    wm /= wsum
                    wb /= wsum
                    wg /= wsum
                weight_combos.append((wm, wb, wg))

    print(f"\nRobustness analysis: varying weights +/-{weight_var*100:.0f}%  "
          f"({len(weight_combos)} weight combos)")

    for r in top_results:
        scores = np.array([
            max(r['mem_score'], 1e-6)**wm * max(r['bal_score'], 1e-6)**wb * max(r['eff_score'], 1e-6)**wg
            for wm, wb, wg in weight_combos
        ])
        r['score_mean'] = float(np.mean(scores))
        r['score_std']  = float(np.std(scores))
        # rank range: what is the worst and best rank this config achieves?
        r['score_min']  = float(np.min(scores))
        r['score_max']  = float(np.max(scores))

    # also compute rank stability: for each weight combo, rank all top configs
    # and record each config's rank distribution
    n_top = len(top_results)
    rank_matrix = np.zeros((len(weight_combos), n_top), dtype=int)

    for wi, (wm, wb, wg) in enumerate(weight_combos):
        combo_scores = [
            max(r['mem_score'], 1e-6)**wm * max(r['bal_score'], 1e-6)**wb * max(r['eff_score'], 1e-6)**wg
            for r in top_results
        ]
        order = np.argsort(combo_scores)[::-1]
        for rank_pos, cfg_idx in enumerate(order):
            rank_matrix[wi, cfg_idx] = rank_pos + 1

    for ci, r in enumerate(top_results):
        ranks = rank_matrix[:, ci]
        r['rank_mean'] = float(np.mean(ranks))
        r['rank_min']  = int(np.min(ranks))
        r['rank_max']  = int(np.max(ranks))

    # print table
    hdr = (
        f"{'#':>3s} | {'Nx':>5s} {'Ny':>4s} {'Nz':>4s} | "
        f"{'NxM':>5s} {'NyM':>5s} {'NzM':>5s} | "
        f"{'Score':>6s} {'mean':>6s} {'std':>6s} {'min':>6s} {'max':>6s} | "
        f"{'Rank':>4s} {'Rmin':>4s} {'Rmax':>4s}"
    )
    sep = "-" * len(hdr)

    print(sep)
    print(hdr)
    print(sep)

    for i, r in enumerate(top_results, 1):
        print(
            f"{i:3d} | "
            f"{r['Nx']:5d} {r['Ny']:4d} {r['Nz']:4d} | "
            f"{r['NxMesh']:5d} {r['NyMesh']:5d} {r['NzMesh']:5d} | "
            f"{r['combined']:6.4f} {r['score_mean']:6.4f} {r['score_std']:6.4f} "
            f"{r['score_min']:6.4f} {r['score_max']:6.4f} | "
            f"{r['rank_mean']:4.1f} {r['rank_min']:4d} {r['rank_max']:4d}"
        )

    print(sep)


# ============================================================================
# Section 10: Main
# ============================================================================

def main():
    p = argparse.ArgumentParser(description="Scan tile/mesh configs for PIC shock simulations")

    p.add_argument('--Lx', type=int, default=65536, help='total x cells (default: 65536)')
    p.add_argument('--Ly', type=int, default=128,   help='total y cells (default: 128)')
    p.add_argument('--Lz', type=int, default=128,   help='total z cells (default: 128)')
    p.add_argument('--ppc', type=int, default=2,     help='particles per cell per species (default: 2)')
    p.add_argument('--n_species', type=int, default=2, help='number of particle species (default: 2)')
    p.add_argument('--sigma', type=float, default=3.0, help='magnetization sigma (default: 3.0)')
    p.add_argument('--n_gpus', type=int, default=8,  help='number of GPUs/MPI ranks (default: 8)')
    p.add_argument('--gpu_mem_gb', type=float, default=40.0, help='GPU memory limit in GB (default: 40)')
    p.add_argument('--min_mesh', type=int, default=8, help='minimum NxMesh/NyMesh/NzMesh (default: 8)')
    p.add_argument('--w_mem', type=float, default=0.3, help='memory score weight (default: 0.3)')
    p.add_argument('--w_bal', type=float, default=0.4, help='balance score weight (default: 0.4)')
    p.add_argument('--w_gpu', type=float, default=0.3, help='GPU efficiency weight (default: 0.3)')
    p.add_argument('--top', type=int, default=20,    help='number of top results to show (default: 20)')
    p.add_argument('--weight_var', type=float, default=0.1, help='weight variation fraction for robustness (default: 0.1 = 10%%)')
    p.add_argument('--plot', action='store_true',    help='generate tradeoff plots')

    args = p.parse_args()

    results = scan_configs(args)

    if results:
        print_table(results, top=args.top)
        robustness_analysis(results, args.w_mem, args.w_bal, args.w_gpu,
                            args.top, args.weight_var)

        if args.plot:
            plot_results(results, args.gpu_mem_gb)
    else:
        print("No valid configurations found. Try relaxing --gpu_mem_gb or --min_mesh.")


if __name__ == '__main__':
    main()
