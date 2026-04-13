"""
Pure-Python Hilbert curve implementation.

Generalized Hilbert space-filling curves for 2D and 3D grids with
arbitrary power-of-two dimensions (2^m0 x 2^m1 [x 2^m2]).

Based on Chris Hamilton, Technical Report CS-2006-07,
Faculty of Computer Science, Dalhousie University, Halifax.
"""

import numpy as np

# ---- bit-manipulation primitives -------------------------------------------

def bit(i, k):
    return (i >> k) & 1

def setbit(i, k, value):
    return (i & ~(1 << k)) | (value << k)

def rotl(value, shift, dim):
    return ((value << shift) | (value >> (dim - shift))) & ((1 << dim) - 1)

def rotr(value, shift, dim):
    return ((value >> shift) | (value << (dim - shift))) & ((1 << dim) - 1)

def gc(i):
    return i ^ (i >> 1)

def gcinv(g):
    i = g
    j = 1
    while (1 << j) <= g:
        i = i ^ (g >> j)
        j += 1
    return i

def tsb(i):
    k = 0
    while (i & 1) != 0:
        i = i >> 1
        k += 1
    return k

def direction(i, dim):
    if i == 0:
        return 0
    if (i & 1) != 0:
        return tsb(i) % dim
    return tsb(i - 1) % dim

def entry(i):
    if i == 0:
        return 0
    return gc(2 * ((i - 1) // 2))

def ted(e, d, b, dim):
    return rotr(b ^ e, d + 1, dim)

def tedinv(e, d, b, dim):
    return rotl(b, d + 1, dim) ^ e


# ---- lookup tables (computed once at import) -------------------------------

def _build_luts(dim):
    """Build numpy lookup tables for all helper functions at a given dimension."""
    n = 1 << dim  # 4 for 2D, 8 for 3D

    gcinv_lut     = np.array([gcinv(i)          for i in range(n)], dtype=np.int64)
    entry_lut     = np.array([entry(i)          for i in range(n)], dtype=np.int64)
    direction_lut = np.array([direction(i, dim) for i in range(n)], dtype=np.int64)

    rotl_lut = np.empty((n, dim), dtype=np.int64)
    for v in range(n):
        for s in range(dim):
            rotl_lut[v, s] = rotl(v, s, dim)

    ted_lut = np.empty((n, dim, n), dtype=np.int64)
    for e in range(n):
        for d in range(dim):
            for b in range(n):
                ted_lut[e, d, b] = ted(e, d, b, dim)

    return gcinv_lut, entry_lut, direction_lut, rotl_lut, ted_lut

_gcinv_3d, _entry_3d, _direction_3d, _rotl_3d, _ted_3d = _build_luts(3)
_gcinv_2d, _entry_2d, _direction_2d, _rotl_2d, _ted_2d = _build_luts(2)


# ---- LUT-based forward Hilbert index (2D) ----------------------------------

def _hilbertindex_2d(m, x, y, einit, dinit):
    """2D cubic Hilbert index. x, y are 1-D int64 arrays.

    Returns (h, e, d) as arrays.
    """
    N = x.shape[0]
    h = np.zeros(N, dtype=np.int64)
    e = np.full(N, einit, dtype=np.int64) if np.isscalar(einit) else einit.copy()
    d = np.full(N, dinit, dtype=np.int64) if np.isscalar(dinit) else dinit.copy()

    for i in range(m - 1, -1, -1):
        l = ((y >> i) & 1) * 2 + ((x >> i) & 1)
        l = _ted_2d[e, d, l]
        w = _gcinv_2d[l]
        e = e ^ _rotl_2d[_entry_2d[w], (d + 1) % 2]
        d = (d + _direction_2d[w] + 1) % 2
        h = (h << 2) | w

    return h, e, d


# ---- LUT-based forward Hilbert index (3D) ----------------------------------

def _hilbertindex_3d(m, x, y, z, einit, dinit):
    """3D cubic Hilbert index. x, y, z are 1-D int64 arrays."""
    N = x.shape[0]
    h = np.zeros(N, dtype=np.int64)
    e = np.full(N, einit, dtype=np.int64) if np.isscalar(einit) else einit.copy()
    d = np.full(N, dinit, dtype=np.int64) if np.isscalar(dinit) else dinit.copy()

    for i in range(m - 1, -1, -1):
        l = ((z >> i) & 1) * 4 + ((y >> i) & 1) * 2 + ((x >> i) & 1)
        l = _ted_3d[e, d, l]
        w = _gcinv_3d[l]
        e = e ^ _rotl_3d[_entry_3d[w], (d + 1) % 3]
        d = (d + _direction_3d[w] + 1) % 3
        h = (h << 3) | w

    return h


# ---- general (non-cubic) forward Hilbert index (2D) ------------------------

def _generalhilbertindex_2d(m0, m1, x, y):
    """General 2D Hilbert index; x, y are 1-D int64 arrays.

    Returns (h, einit, dinit) as arrays.
    """
    dinit = 0
    einit = 0

    if m0 >= m1:
        mmin, mmax = m1, m0
        use_x = True
    else:
        mmin, mmax = m0, m1
        use_x = False
        dinit = 1

    localx = x.copy()
    localy = y.copy()
    target = localx if use_x else localy

    h = np.zeros_like(x, dtype=np.int64)
    for i in range(mmax - 1, mmin - 1, -1):
        l = (target >> i) & 1
        h += l * (1 << (i + mmin))
        target -= l * (1 << i)

    if use_x:
        localx = target
    else:
        localy = target

    if mmin > 0:
        hsub, einit_arr, dinit_arr = _hilbertindex_2d(mmin, localx, localy, einit, dinit)
        h += hsub
        return h, einit_arr, dinit_arr
    else:
        N = x.shape[0]
        return h, np.full(N, einit, dtype=np.int64), np.full(N, dinit, dtype=np.int64)


# ---- general (non-cubic) forward Hilbert index (3D) ------------------------

def _generalhilbertindex_3d(m0, m1, m2, x, y, z):
    """General 3D Hilbert index. x, y, z are 1-D int64 arrays."""
    mi = [m0, m1, m2]

    if mi[0] >= mi[1] and mi[0] >= mi[2]:
        dimmax = 0
    elif mi[1] > mi[0] and mi[1] >= mi[2]:
        dimmax = 1
    else:
        dimmax = 2

    if mi[(dimmax + 1) % 3] >= mi[(dimmax + 2) % 3]:
        dimmed = (dimmax + 1) % 3
        dimmin = (dimmax + 2) % 3
    else:
        dimmed = (dimmax + 2) % 3
        dimmin = (dimmax + 1) % 3

    localp = [x, y, z]

    tempp_max = localp[dimmax] >> mi[dimmin]
    tempp_med = localp[dimmed] >> mi[dimmin]

    h2d, einit_arr, dinit_arr = _generalhilbertindex_2d(
        mi[dimmax] - mi[dimmin],
        mi[dimmed] - mi[dimmin],
        tempp_max, tempp_med
    )
    h = h2d * np.int64(1 << (3 * mi[dimmin]))

    if mi[dimmin] > 0:
        mask = np.int64((1 << mi[dimmin]) - 1)
        tp_max = localp[dimmax] & mask
        tp_med = localp[dimmed] & mask
        tp_min = localp[dimmin] & mask
        h += _hilbertindex_3d(mi[dimmin], tp_max, tp_med, tp_min, einit_arr, dinit_arr)

    return h


# ---- cubic Hilbert index inverse (2D) -------------------------------------

def _hilbertindexinv_2d(m, h, einit, dinit):
    """Inverse 2D cubic Hilbert index. Returns (x, y)."""
    e, d = einit, dinit
    x, y = 0, 0
    for i in range(m - 1, -1, -1):
        w = ((bit(h, 2*i + 1)) << 1) + bit(h, 2*i)
        l = gc(w)
        l = tedinv(e, d, l, 2)
        x = setbit(x, i, bit(l, 0))
        y = setbit(y, i, bit(l, 1))
        e = e ^ rotl(entry(w), d + 1, 2)
        d = (d + direction(w, 2) + 1) % 2
    return x, y


# ---- cubic Hilbert index inverse (3D) -------------------------------------

def _hilbertindexinv_3d(m, h, einit, dinit):
    """Inverse 3D cubic Hilbert index. Returns (x, y, z)."""
    e, d = einit, dinit
    x, y, z = 0, 0, 0
    for i in range(m - 1, -1, -1):
        w = ((bit(h, 3*i + 2)) << 2) + ((bit(h, 3*i + 1)) << 1) + bit(h, 3*i)
        l = gc(w)
        l = tedinv(e, d, l, 3)
        x = setbit(x, i, bit(l, 0))
        y = setbit(y, i, bit(l, 1))
        z = setbit(z, i, bit(l, 2))
        e = e ^ rotl(entry(w), d + 1, 3)
        d = (d + direction(w, 3) + 1) % 3
    return x, y, z


# ---- general (non-cubic) Hilbert index inverse (2D) ------------------------

def _generalhilbertindexinv_2d(m0, m1, h):
    """Inverse general 2D Hilbert index. Returns (x, y)."""
    einit = 0
    dinit = 0
    shift = 0
    localh = h

    if m0 >= m1:
        mmin, mmax = m1, m0
        shift_x = True
        dinit = 0
    else:
        mmin, mmax = m0, m1
        shift_x = False
        dinit = 1

    for i in range(mmax + mmin - 1, mmin + mmin - 1, -1):
        l = bit(localh, i)
        shift += l * (1 << (i - mmin))
        localh -= l * (1 << i)

    x, y = _hilbertindexinv_2d(mmin, localh, einit, dinit)

    if shift_x:
        x += shift
    else:
        y += shift

    return x, y


# ---- general (non-cubic) Hilbert index inverse (3D) ------------------------

def _generalhilbertindexinv_3d(m0, m1, m2, h):
    """Inverse general 3D Hilbert index. Returns (x, y, z)."""
    mi = [m0, m1, m2]

    if (mi[0] >= mi[1]) and (mi[0] >= mi[2]):
        dimmax = 0
    elif (mi[1] > mi[0]) and (mi[1] >= mi[2]):
        dimmax = 1
    else:
        dimmax = 2

    if mi[(dimmax + 1) % 3] >= mi[(dimmax + 2) % 3]:
        dimmed = (dimmax + 1) % 3
        dimmin = (dimmax + 2) % 3
    else:
        dimmed = (dimmax + 2) % 3
        dimmin = (dimmax + 1) % 3

    # locate sub-hypercube via 2D inversion on the reduced domain
    localh = h >> (mi[dimmin] * 3)
    coords = [0, 0, 0]
    coords[dimmax], coords[dimmed] = _generalhilbertindexinv_2d(
        mi[dimmax] - mi[dimmin],
        mi[dimmed] - mi[dimmin],
        localh
    )

    # recover e, d by re-running the 2D forward index
    xa = np.array([coords[dimmax]], dtype=np.int64)
    ya = np.array([coords[dimmed]], dtype=np.int64)
    _, ea, da = _generalhilbertindex_2d(
        mi[dimmax] - mi[dimmin],
        mi[dimmed] - mi[dimmin],
        xa, ya
    )
    e, d = int(ea[0]), int(da[0])

    # transform to global frame
    coords[dimmax] *= (1 << mi[dimmin])
    coords[dimmed] *= (1 << mi[dimmin])
    coords[dimmin] = 0

    # cubic inversion in the local sub-hypercube
    localh = h & ((1 << (mi[dimmin] * 3)) - 1)
    tx, ty, tz = _hilbertindexinv_3d(mi[dimmin], localh, e, d)

    coords[dimmax] += tx
    coords[dimmed] += ty
    coords[dimmin] += tz

    return coords[0], coords[1], coords[2]


# ---- public API ------------------------------------------------------------

class Hilbert2D:
    """Generalized 2D Hilbert curve for a 2^m0 x 2^m1 grid."""

    def __init__(self, m0, m1):
        self.m0 = int(m0)
        self.m1 = int(m1)

    def hindex(self, x, y):
        xa = np.array([x], dtype=np.int64)
        ya = np.array([y], dtype=np.int64)
        h, _, _ = _generalhilbertindex_2d(self.m0, self.m1, xa, ya)
        return int(h[0])

    def inv(self, h):
        return _generalhilbertindexinv_2d(self.m0, self.m1, h)


class Hilbert3D:
    """Generalized 3D Hilbert curve for a 2^m0 x 2^m1 x 2^m2 grid."""

    def __init__(self, m0, m1, m2):
        self.m0 = int(m0)
        self.m1 = int(m1)
        self.m2 = int(m2)

    def hindex(self, x, y, z):
        xa = np.array([x], dtype=np.int64)
        ya = np.array([y], dtype=np.int64)
        za = np.array([z], dtype=np.int64)
        h = _generalhilbertindex_3d(self.m0, self.m1, self.m2, xa, ya, za)
        return int(h[0])

    def hindex_vec(self, x, y, z):
        """Vectorized: x, y, z are 1-D int64 numpy arrays."""
        return _generalhilbertindex_3d(self.m0, self.m1, self.m2, x, y, z)

    def inv(self, h):
        return _generalhilbertindexinv_3d(self.m0, self.m1, self.m2, h)
