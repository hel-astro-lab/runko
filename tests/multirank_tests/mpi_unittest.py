"""
Nano unittest library for MPI related tests.
"""

from mpi4py import MPI
from dataclasses import dataclass

_comm = MPI.COMM_WORLD
_comm_size = _comm.Get_size()
_rank = _comm.Get_rank()


def _gather_error_msgs(msg: str):
    msgs = _comm.gather(msg, root=0)

    if _rank == 0:
        concatted_msg = "\n"
        passes = True
        for msg in msgs:
            passes = passes and len(msg) == 0
            concatted_msg += msg

        if not passes:
            raise RuntimeError(concatted_msg)


def assertEqual(a, b) -> None:
    """
    Asserts a == b on every rank.

    Has to be called on every rank.
    """
    ok = a == b

    err_msg = ""
    if not ok:
        err_msg = f"rank {_rank}: {a} != {b}\n"

    _gather_error_msgs(err_msg)

def assertNotEqual(a, b) -> None:
    """
    Asserts a != b on every rank.

    Has to be called on every rank.
    """
    ok = a != b

    err_msg = ""
    if not ok:
        err_msg = f"rank {_rank}: {a} == {b}\n"

    _gather_error_msgs(err_msg)


@dataclass
class DeferredAssertResult:
    """
    Represents result of deferred assertion result.
    This means that functions returning objects of this class
    do not yet assert on failed test, but store the result in the object.

    These can be checked collectively with assertDeferredResults,
    which will sync accross ranks and assert that they are correct.
    """

    error_msg: str = ""

    def fail(self) -> bool:
        return 0 != len(self.error_msg)


def assertEqualDeferred(a, b) -> DeferredAssertResult:
    """
    Similar to assertEqual but,
    the checking is deferred to assertDeferredResults.
    """

    if a == b:
        return DeferredAssertResult("")
    else:
        return DeferredAssertResult(f"rank {_rank}: {a} != {b}\n")


def assertDeferredResults(results: list[DeferredAssertResult]) -> None:
    """
    Assert that all given results pass on every rank,
    i.e. has to be called on every rank.
    """

    concatted_err_msg = ""

    for result in results:
        concatted_err_msg += result.error_msg

    _gather_error_msgs(concatted_err_msg)
