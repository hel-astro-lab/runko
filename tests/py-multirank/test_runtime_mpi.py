# Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

import mpi_unittest as mpiut
import runko


if __name__ == "__main__":
    runtime = runko.actions.RuntimeInstance()
    rank, size = runtime.mpi_rank(), runtime.mpi_size()

    mpiut.assertEqual(True, 2 <= size)
    mpiut.assertEqual(rank, runko.default_runtime().mpi_rank())
    mpiut.assertEqual(runko.on_main_rank(), rank == 0)

    runtime.mpi_barrier()

    for root in range(size):
        payload = {"root": root, "data": list(range(root + 3))} if rank == root else None
        received = runtime.mpi_bcast(payload, root=root)
        mpiut.assertEqual(received, {"root": root, "data": list(range(root + 3))})


    for root in range(size):
        # Different sized messages from each rank.
        gathered = runtime.mpi_gather("x" * rank, root=root)
        expected = ["x" * i for i in range(size)] if rank == root else None
        mpiut.assertEqual(gathered, expected)


    ranks = runtime.mpi_gather(rank)
    if rank == 0:
        mpiut.assertEqual(ranks, list(range(size)))
    else:
        mpiut.assertEqual(ranks, None)
