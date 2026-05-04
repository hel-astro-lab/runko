# Copyright 2025 - 2026, Joonas Nättilä and the hel-astro-lab contributors
# SPDX-License-Identifier: GPL-3.0-or-later

import mpi_unittest
import random

if __name__ == "__main__":
    mpi_unittest.assertEqual(45, 45)

    results = []
    for i in range(random.randint(0, 10)):
        results.append(mpi_unittest.assertEqualDeferred(int, int))

    mpi_unittest.assertDeferredResults(results)
