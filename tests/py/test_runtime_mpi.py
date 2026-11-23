# Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

import sys
import unittest

import runko
import runko.runtime


class runtime_mpi(unittest.TestCase):
    def test_import_does_not_initialize_runtime(self):
        # Checked in __main__ before anything creates a runtime.
        self.assertFalse(_runtime_existed_after_import)
        self.assertNotIn("mpi4py", sys.modules)


    def test_rank_and_size(self):
        r = runko.actions.RuntimeInstance()
        self.assertEqual(r.mpi_rank(), 0)
        self.assertEqual(r.mpi_size(), 1)
        self.assertTrue(runko.on_main_rank())


    def test_default_runtime_is_cached(self):
        self.assertIs(runko.default_runtime(), runko.default_runtime())


    def test_barrier(self):
        runko.default_runtime().mpi_barrier()


    def test_bcast(self):
        r = runko.default_runtime()
        for obj in (None, 1, "foo", "", {"a": [1, 2.5, (3, "x")]}, b"\x00\x01"):
            self.assertEqual(r.mpi_bcast(obj), obj)
            self.assertEqual(r.mpi_bcast(obj, root=0), obj)


    def test_gather(self):
        r = runko.default_runtime()
        self.assertEqual(r.mpi_gather("foo"), ["foo"])
        self.assertEqual(r.mpi_gather({"x": 1}, root=0), [{"x": 1}])


    def test_invalid_root(self):
        r = runko.default_runtime()
        with self.assertRaises(ValueError):
            r.mpi_bcast(1, root=1)
        with self.assertRaises(ValueError):
            r.mpi_gather(1, root=-1)




_runtime_existed_after_import = runko.runtime._default_runtime is not None

if __name__ == "__main__":
    unittest.main()
