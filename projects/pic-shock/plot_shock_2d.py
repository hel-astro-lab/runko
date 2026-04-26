import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import argparse
import math
import os

import runko
from runko.auto_outdir import resolve_outdir
from runko.mpiio_reader import read_field_snapshot
from runko.mpiio_spectra_reader import (
    read_spectra_header,
    read_spectra_snapshot,
    u_bin_edges,
    beta_bin_edges,
)


# -- variable configuration ---------------------------------------------------

VAR_CONFIG = {
    'rho': {
        'cmap': 'viridis',
        'vmin': 0.0,
        'vmax': 5.0,
        'clabel': r'$n/n_0$',
        'symmetric': False,
        'norm_key': 'rho',
    },
    'ex': {
        'cmap': 'PiYG',
        'vmin': None,
        'vmax': None,
        'clabel': r'$E_x/B_0$',
        'symmetric': True,
        'norm_key': 'ex',
    },
    'ey': {
        'cmap': 'PiYG',
        'vmin': None,
        'vmax': None,
        'clabel': r'$E_y/B_0$',
        'symmetric': True,
        'norm_key': 'ey',
    },
    'ez': {
        'cmap': 'PiYG',
        'vmin': None,
        'vmax': None,
        'clabel': r'$E_z/B_0$',
        'symmetric': True,
        'norm_key': 'ez',
    },
    'bx': {
        'cmap': 'BrBG',
        'vmin': None,
        'vmax': None,
        'clabel': r'$B_x/B_0$',
        'symmetric': True,
        'norm_key': 'bx',
    },
    'by': {
        'cmap': 'BrBG',
        'vmin': None,
        'vmax': None,
        'clabel': r'$B_y/B_0$',
        'symmetric': True,
        'norm_key': 'by',
    },
    'bz': {
        'cmap': 'BrBG',
        'vmin': None,
        'vmax': None,
        'clabel': r'$B_z/B_0$',
        'symmetric': True,
        'norm_key': 'bz',
    },
    'jx': {
        'cmap': 'RdBu',
        'vmin': None,
        'vmax': None,
        'clabel': r'$J_x/(n_0 q_0 c)$',
        'symmetric': True,
        'norm_key': 'jx',
    },
    'jy': {
        'cmap': 'RdBu',
        'vmin': None,
        'vmax': None,
        'clabel': r'$J_y/(n_0 q_0 c)$',
        'symmetric': True,
        'norm_key': 'jy',
    },
    'jz': {
        'cmap': 'RdBu',
        'vmin': None,
        'vmax': None,
        'clabel': r'$J_z/(n_0 q_0 c)$',
        'symmetric': True,
        'norm_key': 'jz',
    },
    's0_u': {
        'cmap': 'plasma_r',
        'vmin': -2,
        'vmax': 5,
        'clabel': r'$\log_{10} f_e(u)$',
        'ylabel': r'$\log_{10} u$',
        'symmetric': False,
    },
    's0_bx': {
        'cmap': 'plasma_r',
        'vmin': -2,
        'vmax': 5,
        'clabel': r'$\log_{10} f_e(\beta_x)$',
        'ylabel': r'$\beta_x$',
        'symmetric': False,
    },
    's0_by': {
        'cmap': 'plasma_r',
        'vmin': -2,
        'vmax': 5,
        'clabel': r'$\log_{10} f_e(\beta_y)$',
        'ylabel': r'$\beta_y$',
        'symmetric': False,
    },
    's0_bz': {
        'cmap': 'plasma_r',
        'vmin': -2,
        'vmax': 5,
        'clabel': r'$\log_{10} f_e(\beta_z)$',
        'ylabel': r'$\beta_z$',
        'symmetric': False,
    },
    's1_u': {
        'cmap': 'plasma_r',
        'vmin': -2,
        'vmax': 5,
        'clabel': r'$\log_{10} f_i(u)$',
        'ylabel': r'$\log_{10} u$',
        'symmetric': False,
    },
    's1_bx': {
        'cmap': 'plasma_r',
        'vmin': -2,
        'vmax': 5,
        'clabel': r'$\log_{10} f_i(\beta_x)$',
        'ylabel': r'$\beta_x$',
        'symmetric': False,
    },
    's1_by': {
        'cmap': 'plasma_r',
        'vmin': -2,
        'vmax': 5,
        'clabel': r'$\log_{10} f_i(\beta_y)$',
        'ylabel': r'$\beta_y$',
        'symmetric': False,
    },
    's1_bz': {
        'cmap': 'plasma_r',
        'vmin': -2,
        'vmax': 5,
        'clabel': r'$\log_{10} f_i(\beta_z)$',
        'ylabel': r'$\beta_z$',
        'symmetric': False,
    },
}


# -- physics helpers (adapted from plot3d_pyvista_shock.py) --------------------

def add_shock_derived(conf):
    """Compute derived charge/mass/B-field quantities on a v5 Configuration."""
    oppc = 2 * conf.ppc
    omp  = conf.cfl / conf.c_omp
    conf.qe = -(omp**2 * conf.upstream_gamma) / (0.5 * oppc * (1.0 + conf.m0 / conf.m1))
    conf.me = conf.m0 * abs(conf.qe)
    conf.binit = math.sqrt(
        conf.upstream_gamma * oppc * 0.5 * conf.cfl**2
        * conf.me * (1.0 + conf.me / (conf.m1 * abs(conf.qe)))
        * conf.sigma
    )


def get_normalization(var, conf):
    """Return the normalization factor for a given field variable.

    The MPI-IO field dump volume-sums n0,n1,jx,jy,jz over the stride^3 fine
    cells that map to one coarse output cell, so number/current densities
    pick up an extra stride^3 factor. E,B are subsampled (no stride factor).
    """
    stride3 = conf.stride**3
    n0_coarse = 2.0 * conf.ppc * stride3   # upstream coarse-cell density
    qe = abs(conf.qe)

    if var == 'rho':
        return n0_coarse
    if var in ('jx', 'jy', 'jz'):
        # J ~ n_0 q_0 c (with c -> cfl in code units), summed over stride^3
        return n0_coarse * qe * conf.cfl
    if var in ('bx', 'by', 'bz', 'ex', 'ey', 'ez'):
        return conf.binit
    return 1.0


# -- shock tracking ------------------------------------------------------------

def find_shock_index(rho_1d, threshold):
    """Locate shock front by walking from the right (high x) leftward.

    Returns the largest x-index where rho_1d > threshold; 0 if no cell
    exceeds the threshold (shock not yet formed).
    """
    above = np.where(rho_1d > threshold)[0]
    if len(above) == 0:
        return 0
    return int(above[-1])


# -- plotting helpers ----------------------------------------------------------

def winsorize_clim(data, lo=0.005, hi=0.005):
    """Return (vmin, vmax) from winsorized quantiles."""
    return np.quantile(data, [lo, 1.0 - hi])


def add_colorbar(fig, ax, im, clabel):
    """Add a thin vertical colorbar to the right of ax."""
    pos = ax.get_position()
    pad = 0.008
    cax = fig.add_axes([pos.x1 + pad, pos.y0, 0.008, pos.height])
    cb = fig.colorbar(im, cax=cax, orientation='vertical', ticklocation='right')
    cax.tick_params(labelsize=5)
    cax.set_ylabel(clabel, fontsize=6, rotation=270, labelpad=8)
    return cb


def plot_density_profile(ax, rho_1d, x_coords, cfg):
    """Plot 1D y,z-averaged density profile."""
    ax.plot(x_coords, rho_1d, color='C0', linewidth=0.8)
    ax.set_ylabel(cfg['clabel'])
    ax.set_ylim(cfg['vmin'], cfg['vmax'])
    ax.minorticks_on()


def plot_field_panel(ax, fig, data_2d, xmin, xmax, zmin, zmax, cfg):
    """Plot a 2D field slice with imshow and colorbar."""
    vmin, vmax = cfg['vmin'], cfg['vmax']

    if vmin is None or vmax is None:
        wmin, wmax = winsorize_clim(data_2d)
        vmin = wmin if vmin is None else vmin
        vmax = wmax if vmax is None else vmax

    if cfg['symmetric']:
        vlim = max(abs(vmin), abs(vmax))
        vmin, vmax = -vlim, vlim

    im = ax.imshow(
        data_2d,
        origin='lower',
        aspect='auto',
        extent=[xmin, xmax, zmin, zmax],
        cmap=cfg['cmap'],
        vmin=vmin,
        vmax=vmax,
        interpolation='nearest',
    )
    ax.minorticks_on()
    add_colorbar(fig, ax, im, cfg['clabel'])
    return im


def plot_spectra_panel(ax, fig, spec_2d, x_edges, bin_edges, cfg,
                       log_y=False, dark=False):
    """Plot a 2D spectrum (x vs momentum bin) with pcolormesh and colorbar."""
    vmin, vmax = cfg['vmin'], cfg['vmax']

    if log_y:
        bin_plot = np.log10(bin_edges)
    else:
        bin_plot = bin_edges

    # mask bins with zero counts (log10 of zero → -inf)
    masked = np.ma.masked_where(spec_2d < vmin, spec_2d)

    bg = 'black' if dark else 'white'
    cmap = matplotlib.colormaps[cfg['cmap']].copy()
    cmap.set_bad(color=bg)
    cmap.set_under(color=bg)
    ax.set_facecolor(bg)

    im = ax.pcolormesh(
        x_edges,
        bin_plot,
        masked.T,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        shading='flat',
    )
    ax.minorticks_on()

    if not log_y:
        ax.axhline(y=0.0, linestyle='dotted', alpha=0.4, color='grey', linewidth=0.5)

    add_colorbar(fig, ax, im, cfg['clabel'])
    return im


# -- main ----------------------------------------------------------------------

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="2D multi-panel shock visualization")
    parser.add_argument("--conf", required=True, help="Path to .ini config file")
    parser.add_argument("--lap", type=int, required=True, help="Lap number")
    parser.add_argument("--dark", action="store_true", help="Use dark background style")
    parser.add_argument("--track_shock", action="store_true",
                        help="Crop the figure to a window locked on the shock front")
    parser.add_argument("--track_shock_threshold", type=float, default=1.5,
                        help="density threshold n/n_0 used to locate the shock front")
    parser.add_argument("--track_shock_win_neg", type=float, default=30.0,
                        help="skin depths to plot behind (downstream of) the shock")
    parser.add_argument("--track_shock_win_pos", type=float, default=30.0,
                        help="skin depths to plot ahead (upstream) of the shock")
    args_cli = parser.parse_args()

    # -- config ----------------------------------------------------------------
    conf = runko.Configuration(args_cli.conf)
    add_shock_derived(conf)
    conf.outdir = resolve_outdir(conf)
    lap = args_cli.lap

    # output cell width in skin depths: each dump cell spans `stride` fine cells.
    dx = conf.stride / conf.c_omp

    # -- style -----------------------------------------------------------------
    if args_cli.dark: plt.style.use('dark_background')

    plt.rc('font', family='serif')
    plt.rc('text', usetex=False)
    plt.rc('xtick', top=True, direction='out', labelsize=7)
    plt.rc('ytick', right=True, direction='out', labelsize=7)
    plt.rc('axes', labelsize=8)

    # -- read field data -------------------------------------------------------
    flds_path = os.path.join(conf.outdir, f"flds_{lap}.bin")
    print(f"reading {flds_path}")
    fields = read_field_snapshot(flds_path)

    nz, ny, nx = fields['ex'].shape
    print(f"field shape (nz, ny, nx) = ({nz}, {ny}, {nx})")

    # coordinate arrays (skin depths)
    x_coords = np.arange(nx) * dx
    z_coords = np.arange(nz) * dx
    xmin, xmax = x_coords[0], x_coords[-1]
    zmin, zmax = z_coords[0], z_coords[-1]

    # total density, normalized
    norm_rho = get_normalization('rho', conf)
    rho = (fields['n0'] + fields['n1']) / norm_rho

    # 1D density: average over y,z → shape (nx,)
    rho_1d = rho.mean(axis=(0, 1))

    # 2D midplane slices (y = ny//2) → shape (nz, nx)
    ymid = ny // 2
    rho_2d = rho[:, ymid, :]

    field_vars = ['ex', 'ey', 'ez', 'bx', 'by', 'bz', 'jx', 'jy', 'jz']
    field_slices = {}
    for v in field_vars:
        norm = get_normalization(v, conf)
        field_slices[v] = fields[v][:, ymid, :] / norm

    # -- read spectra data (optional) -----------------------------------------
    spec_panel_order = ['s0_u', 's0_bx', 's0_by', 's0_bz',
                        's1_u', 's1_bx', 's1_by', 's1_bz']
    spec_path = os.path.join(conf.outdir, f"pspectra_{lap}.bin")
    spec_data = {}
    spec_x_edges = None
    u_edges = b_edges = None
    if os.path.exists(spec_path):
        print(f"reading {spec_path}")
        spec_hdr = read_spectra_header(spec_path)
        spectra = read_spectra_snapshot(spec_path)
        u_edges = u_bin_edges(spec_hdr)
        b_edges = beta_bin_edges(spec_hdr)
        spec_nx = spec_hdr['nx']
        spec_x_edges = np.arange(spec_nx + 1) * dx
        for sv in spec_panel_order:
            if sv in spectra:
                raw = spectra[sv].sum(axis=(0, 1))  # (nx, nbins)
                spec_data[sv] = np.log10(raw + 1e-30)
    else:
        print(f"no spectra file at {spec_path}, skipping spectra panels")

    # -- optional: lock view on the shock front --------------------------------
    xlabel = r'$x~(c/\omega_p)$'
    if args_cli.track_shock:
        ix_shock = find_shock_index(rho_1d, args_cli.track_shock_threshold)
        x_shock = ix_shock * dx
        if ix_shock == 0:
            print(f"no cell exceeds n/n_0 > {args_cli.track_shock_threshold}; "
                  f"locking shock to x=0")
        else:
            print(f"shock detected at ix={ix_shock}, x={x_shock:.2f} c/omega_p")

        # window edges in skin depths; convert to cell indices and clamp.
        ix_lo = max(0,  int(round((x_shock - args_cli.track_shock_win_neg) / dx)))
        ix_hi = min(nx, int(round((x_shock + args_cli.track_shock_win_pos) / dx)) + 1)
        if ix_hi <= ix_lo:
            ix_lo, ix_hi = 0, max(1, min(nx, int(round(args_cli.track_shock_win_pos / dx)) + 1))

        # slice 1D and 2D field arrays; rebuild shock-relative x coordinates
        rho_1d  = rho_1d[ix_lo:ix_hi]
        rho_2d  = rho_2d[:, ix_lo:ix_hi]
        for v in field_vars:
            field_slices[v] = field_slices[v][:, ix_lo:ix_hi]
        x_coords = np.arange(ix_lo, ix_hi) * dx - x_shock
        xmin, xmax = x_coords[0], x_coords[-1]

        # spectra: shift x edges to shock-relative; slice if dump nx matches.
        if spec_data and spec_x_edges is not None:
            if spec_x_edges.shape[0] - 1 == nx:
                spec_x_edges = spec_x_edges[ix_lo:ix_hi + 1] - x_shock
                for sv in list(spec_data.keys()):
                    spec_data[sv] = spec_data[sv][ix_lo:ix_hi]
            else:
                spec_x_edges = spec_x_edges - x_shock

        xlabel = r'$x - x_\mathrm{sh}~(c/\omega_p)$'

    # -- figure setup ----------------------------------------------------------
    n_field_panels = 2 + len(field_vars)            # 1D rho + 2D rho + E,B,J
    n_spec_panels  = len(spec_data)
    npanels = n_field_panels + n_spec_panels
    fig_height = 1.5 + 1.0 * (npanels - 1)
    fig = plt.figure(1, figsize=(7.0, fig_height))
    height_ratios = [1.5] + [1.0] * (npanels - 1)
    gs = plt.GridSpec(npanels, 1, hspace=0, height_ratios=height_ratios)

    axs = []
    for i in range(npanels):
        axs.append(plt.subplot(gs[i, 0]))
        axs[i].minorticks_on()

    fig.subplots_adjust(left=0.10, right=0.85, bottom=0.025, top=0.995)

    # -- panel 0: 1D density profile ------------------------------------------
    plot_density_profile(axs[0], rho_1d, x_coords, VAR_CONFIG['rho'])

    # -- panel 1: 2D density --------------------------------------------------
    plot_field_panel(axs[1], fig, rho_2d, xmin, xmax, zmin, zmax, VAR_CONFIG['rho'])

    # -- panels 2..2+len(field_vars)-1: E, B, J field slices -------------------
    for i, v in enumerate(field_vars):
        plot_field_panel(axs[2 + i], fig, field_slices[v], xmin, xmax, zmin, zmax, VAR_CONFIG[v])

    # -- spectra panels (only those present) ----------------------------------
    spec_axes_start = n_field_panels
    spec_idx = 0
    for sv in spec_panel_order:
        if sv not in spec_data:
            continue
        is_u = sv.endswith('_u')
        edges = u_edges if is_u else b_edges
        ax = axs[spec_axes_start + spec_idx]
        plot_spectra_panel(ax, fig, spec_data[sv], spec_x_edges, edges,
                           VAR_CONFIG[sv], log_y=is_u, dark=args_cli.dark)
        ax.set_ylabel(VAR_CONFIG[sv]['ylabel'])
        spec_idx += 1

    # -- shared x-axis formatting ----------------------------------------------
    for ax in axs[:-1]:
        ax.tick_params(labelbottom=False)
    axs[-1].set_xlabel(xlabel)

    # set consistent x limits
    for ax in axs:
        ax.set_xlim(xmin, xmax)

    # z-axis label for all 2D field panels (everything between rho-2D and spectra)
    for ax in axs[1:n_field_panels]:
        ax.set_ylabel(r'$z~(c/\omega_p)$')

    # -- save ------------------------------------------------------------------
    slap = str(lap).rjust(7, '0')
    suffix = "_tracked" if args_cli.track_shock else ""
    fname_base = os.path.join(conf.outdir, f"shock_2d_{slap}{suffix}")

    #plt.savefig(fname_base + '.pdf')
    #print(f"saved {fname_base}.pdf")

    plt.savefig(fname_base + '.png', dpi=300)
    print(f"saved {fname_base}.png")
