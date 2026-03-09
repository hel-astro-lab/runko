"""Auto-naming utility for runko output directories."""

import math


def _compact(v):
    """Format a number compactly: drop dots, strip trailing zeros, keep sci notation."""

    # integers pass through
    if isinstance(v, int):
        return str(v)

    # floats that are whole numbers become ints
    if isinstance(v, float) and v == math.floor(v) and abs(v) < 1e15:
        return str(int(v))

    # use scientific notation for very small values that :g would expand
    if abs(v) < 0.01 and v != 0:
        # format in scientific notation, strip trailing zeros from coefficient
        s = f"{v:e}"
        coeff, exp = s.split("e")
        coeff = coeff.rstrip("0").rstrip(".")
        coeff_f = float(coeff)
        exp_i = int(exp)

        if "." not in coeff:
            return f"{int(coeff_f)}e{exp_i:+03d}"

        # shift decimal: 2.5e-03 -> 25e-04
        frac = coeff.split(".")[1]
        shift = len(frac)
        new_coeff = int(round(coeff_f * 10**shift))
        new_exp = exp_i - shift
        return f"{new_coeff}e{new_exp:+03d}"

    # regular decimal: drop the dot, strip trailing zeros
    s = f"{v:g}"
    s = s.replace(".", "")
    return s


def resolve_outdir(config):
    """Resolve the output directory name from a Configuration.

    - Explicit string (not "auto") → returned as-is
    - None / falsy                 → "runko_output"
    - "auto"                       → auto-generated name from config params
    """

    outdir = config.outdir

    if not outdir:
        return "runko_output"

    if outdir != "auto":
        return outdir

    # --- auto mode ---
    parts = []

    # grid dimensions (mandatory for auto)
    gx = config.Nx * config.NxMesh
    gy = config.Ny * config.NyMesh
    gz = config.Nz * config.NzMesh
    parts.append(f"{gx}x{gy}x{gz}")

    # optional tagged params in fixed order
    _tag_params = [
        ("ppc",  "ppc"),
        ("c",    "c_omp"),
        ("s",    "sigma"),
        ("np",   "n_filter_passes"),
        ("cfl",  "cfl"),
        ("t",    "theta"),
        ("gam",  "upstream_gamma"),
    ]

    for tag, attr in _tag_params:
        val = getattr(config, attr, None)
        if val is not None:
            parts.append(f"{tag}{_compact(val)}")

    # B-field projections (all three or nothing)
    bx = config.bx_proj
    by = config.by_proj
    bz = config.bz_proj
    if bx is not None and by is not None and bz is not None:
        parts.append(f"bx{_compact(bx)}by{_compact(by)}bz{_compact(bz)}")

    name = "_".join(parts)

    # prefix/postfix are user-controlled (user includes own underscores)
    if config.prefix is not None:
        name = config.prefix + name
    if config.postfix is not None:
        name = name + config.postfix

    # collapse any accidental double underscores
    while "__" in name:
        name = name.replace("__", "_")

    return name
