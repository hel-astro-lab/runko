import unittest

from runko.auto_outdir import _compact, resolve_outdir
from runko.configuration import Configuration


class TestCompact(unittest.TestCase):
    """Tests for the _compact() number formatter."""

    def test_int_passthrough(self):
        self.assertEqual(_compact(3), "3")
        self.assertEqual(_compact(0), "0")
        self.assertEqual(_compact(100), "100")

    def test_float_to_int(self):
        """Floats with no fractional part become ints."""
        self.assertEqual(_compact(3.0), "3")
        self.assertEqual(_compact(1.0), "1")

    def test_decimal_dot_dropped(self):
        """Decimal point is removed, keeping all significant digits."""
        self.assertEqual(_compact(0.1), "01")
        self.assertEqual(_compact(0.45), "045")

    def test_trailing_zeros_stripped(self):
        self.assertEqual(_compact(1.50), "15")
        self.assertEqual(_compact(0.10), "01")

    def test_scientific_notation(self):
        self.assertEqual(_compact(1e-5), "1e-05")
        self.assertEqual(_compact(2.5e-3), "25e-04")


class TestResolveOutdir(unittest.TestCase):
    """Tests for the resolve_outdir() function."""

    def _make_config(self, **kwargs):
        """Create a Configuration with given attributes."""
        conf = Configuration(None)
        for k, v in kwargs.items():
            setattr(conf, k, v)
        return conf

    def test_explicit_outdir(self):
        """Explicit outdir is returned unchanged."""
        conf = self._make_config(outdir="my_results")
        self.assertEqual(resolve_outdir(conf), "my_results")

    def test_none_defaults(self):
        """None outdir defaults to 'runko_output'."""
        conf = self._make_config()
        self.assertEqual(resolve_outdir(conf), "runko_output")

    def test_auto_all_params(self):
        conf = self._make_config(
            outdir="auto",
            Nx=8, Ny=2, Nz=2,
            NxMesh=32, NyMesh=32, NzMesh=32,
            ppc=2,
            c_omp=3,
            sigma=0.1,
            n_filter_passes=4,
            cfl=0.45,
            theta=1e-5,
            upstream_gamma=3.0,
            bx_proj=0.0, by_proj=0.0, bz_proj=1.0,
        )
        result = resolve_outdir(conf)
        self.assertEqual(result,
                         "256x64x64_ppc2_c3_s01_np4_cfl045_t1e-05_gam3_bx0by0bz1")

    def test_auto_missing_optional(self):
        """Missing optional params are simply skipped."""
        conf = self._make_config(
            outdir="auto",
            Nx=4, Ny=4, Nz=1,
            NxMesh=16, NyMesh=16, NzMesh=16,
            ppc=4,
        )
        result = resolve_outdir(conf)
        self.assertEqual(result, "64x64x16_ppc4")

    def test_auto_prefix_only(self):
        conf = self._make_config(
            outdir="auto",
            prefix="3d_",
            Nx=8, Ny=2, Nz=2,
            NxMesh=32, NyMesh=32, NzMesh=32,
        )
        result = resolve_outdir(conf)
        self.assertEqual(result, "3d_256x64x64")

    def test_auto_postfix_only(self):
        conf = self._make_config(
            outdir="auto",
            postfix="_run7",
            Nx=2, Ny=2, Nz=2,
            NxMesh=10, NyMesh=10, NzMesh=10,
        )
        result = resolve_outdir(conf)
        self.assertEqual(result, "20x20x20_run7")

    def test_double_underscore_collapsed(self):
        """Double underscores from prefix/postfix are collapsed."""
        conf = self._make_config(
            outdir="auto",
            prefix="run__",
            postfix="__end",
            Nx=2, Ny=2, Nz=2,
            NxMesh=10, NyMesh=10, NzMesh=10,
        )
        result = resolve_outdir(conf)
        self.assertNotIn("__", result)
        self.assertEqual(result, "run_20x20x20_end")

    def test_bfield_partial_omitted(self):
        """If only some B-field projections are given, skip the B-field tag."""
        conf = self._make_config(
            outdir="auto",
            Nx=2, Ny=2, Nz=2,
            NxMesh=10, NyMesh=10, NzMesh=10,
            bx_proj=1.0,
        )
        result = resolve_outdir(conf)
        self.assertEqual(result, "20x20x20")

    def test_integration_full_config(self):
        """Full config with prefix, postfix, and all optional params."""
        conf = self._make_config(
            outdir="auto",
            prefix="3d_",
            postfix="_v1_test1",
            Nx=8, Ny=2, Nz=2,
            NxMesh=32, NyMesh=32, NzMesh=32,
            ppc=2,
            c_omp=3,
            sigma=0.1,
            n_filter_passes=4,
            cfl=0.45,
            theta=1e-5,
            upstream_gamma=3.0,
            bx_proj=0.0, by_proj=0.0, bz_proj=1.0,
        )
        result = resolve_outdir(conf)
        self.assertEqual(result,
                         "3d_256x64x64_ppc2_c3_s01_np4_cfl045_t1e-05_gam3_bx0by0bz1_v1_test1")


if __name__ == "__main__":
    unittest.main()
