"""
Nano unittest library for MPI related tests.
"""

from mpi4py import MPI

_comm = MPI.COMM_WORLD
_comm_size = _comm.Get_size()
_rank = _comm.Get_rank()

def assertEqual(a, b):
    """
    Asserts a == b on every rank.
    """

    ok = a == b

    err_msg = ""
    if not ok:
        err_msg = f"rank {_rank}: {a} != {b}\n"

    msgs = _comm.gather(err_msg, root=0)
    if _rank == 0:
        concatted_msg = "\n"
        passes = True
        for msg in msgs:
            passes = passes and len(msg) == 0
            concatted_msg += msg

        if not passes:
            raise RuntimeError(concatted_msg)

