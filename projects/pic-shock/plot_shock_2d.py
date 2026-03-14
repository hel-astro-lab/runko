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
        'clabel': r'$E_x$',
        'symmetric': True,
        'norm_key': 'ex',
    },
    'ey': {
        'cmap': 'PiYG',
        'vmin': None,
        'vmax': None,
        'clabel': r'$E_y$',
        'symmetric': True,
        'norm_key': 'ey',
    },
    'ez': {
        'cmap': 'PiYG',
        'vmin': None,
        'vmax': None,
        'clabel': r'$E_z$',
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
        'clabel': r'$J_x$',
        'symmetric': True,
        'norm_key': 'jx',
    },
    'jy': {
        'cmap': 'RdBu',
        'vmin': None,
        'vmax': None,
        'clabel': r'$J_y$',
        'symmetric': True,
        'norm_key': 'jy',
    },
    'jz': {
        'cmap': 'RdBu',
        'vmin': None,
        'vmax': None,
        'clabel': r'$J_z$',
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
    """Return the normalization factor for a given field variable."""
    n0 = 2.0 * conf.ppc
    qe = abs(conf.qe)
    me_per_qe = abs(conf.me) / qe

    if var == 'rho':
        return n0
    if var in ('jx', 'jy', 'jz'):
        return qe * n0 * conf.cfl**2
    if var in ('bx', 'by', 'bz'):
        return conf.binit
    if var in ('ex', 'ey', 'ez'):
        deltax = 1.0 / conf.c_omp
        return (me_per_qe * conf.cfl**2) / deltax
    return 1.0


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
    args_cli = parser.parse_args()

    # -- config ----------------------------------------------------------------
    conf = runko.Configuration(args_cli.conf)
    add_shock_derived(conf)
    conf.outdir = resolve_outdir(conf)
    lap = args_cli.lap

    dx = 1.0 / conf.c_omp  # cell size in skin depths

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

    # -- read spectra data -----------------------------------------------------
    spec_path = os.path.join(conf.outdir, f"pspectra_{lap}.bin")
    print(f"reading {spec_path}")
    spec_hdr = read_spectra_header(spec_path)
    spectra = read_spectra_snapshot(spec_path)

    # bin edges for pcolormesh
    u_edges  = u_bin_edges(spec_hdr)
    b_edges  = beta_bin_edges(spec_hdr)

    # x edges for spectra (spectra nx may differ from field nx due to stride)
    spec_nx = spec_hdr['nx']
    spec_dx = dx  # same cell-to-skindepth ratio
    spec_x_edges = np.arange(spec_nx + 1) * spec_dx

    # sum over y,z tile axes → (nx, nbins), then log10
    spec_vars = ['s0_u', 's0_bx', 's0_by', 's0_bz',
                 's1_u', 's1_bx', 's1_by', 's1_bz']
    spec_data = {}
    for sv in spec_vars:
        if sv in spectra:
            raw = spectra[sv].sum(axis=(0, 1))  # (nx, nbins)
            spec_data[sv] = np.log10(raw + 1e-30)

    # -- figure setup ----------------------------------------------------------
    npanels = 19
    fig = plt.figure(1, figsize=(7.0, 22.0))
    height_ratios = [1.5] + [1.0] * 18  # 1D density panel is taller
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

    # -- panels 2-10: E, B, J field slices ------------------------------------
    for i, v in enumerate(field_vars):
        plot_field_panel(axs[2 + i], fig, field_slices[v], xmin, xmax, zmin, zmax, VAR_CONFIG[v])

    # -- panels 11-14: species 0 spectra --------------------------------------
    spec_panel_order = ['s0_u', 's0_bx', 's0_by', 's0_bz',
                        's1_u', 's1_bx', 's1_by', 's1_bz']
    for i, sv in enumerate(spec_panel_order):
        if sv not in spec_data:
            continue
        is_u = sv.endswith('_u')
        edges = u_edges if is_u else b_edges
        plot_spectra_panel(axs[11 + i], fig, spec_data[sv], spec_x_edges, edges,
                           VAR_CONFIG[sv], log_y=is_u, dark=args_cli.dark)

        axs[11 + i].set_ylabel(VAR_CONFIG[sv]['ylabel'])

    # -- shared x-axis formatting ----------------------------------------------
    for ax in axs[:-1]:
        ax.tick_params(labelbottom=False)
    axs[-1].set_xlabel(r'$x~(c/\omega_p)$')

    # set consistent x limits
    for ax in axs:
        ax.set_xlim(xmin, xmax)

    # z-axis label for all 2D field panels
    for ax in axs[1:11]:
        ax.set_ylabel(r'$z~(c/\omega_p)$')

    # -- save ------------------------------------------------------------------
    slap = str(lap).rjust(7, '0')
    fname_base = os.path.join(conf.outdir, f"shock_2d_{slap}")

    #plt.savefig(fname_base + '.pdf')
    #print(f"saved {fname_base}.pdf")

    plt.savefig(fname_base + '.png', dpi=300)
    print(f"saved {fname_base}.png")
