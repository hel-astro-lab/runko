"""
Internal logging utilities which are configurable through standard logger module.

Logging in runko is done using named loggers `runko-N`,
where N is the rank where the logging comes from.

By default only main rank will have default logging level changed to INFO
while other ranks have logging level WARNING (python default).
The default handler used in runko is `logging.StreamHandler`
which logging level is set to `logging.NOTSET`.
"""


import logging
from mpi4py import MPI


def on_main_rank() -> bool:
    """
    Checks if the caller is "main" rank.

    There is no other guarantees other that there is only one main rank.
    """

    return MPI.COMM_WORLD.Get_rank() == 0


def runko_logger_name(sublogger=None):
    name = f"runko-{MPI.COMM_WORLD.Get_rank()}"
    if sublogger:
        name += f".{sublogger}"
    return name


def runko_logger(sublogger=None):
    return logging.getLogger(runko_logger_name(sublogger))


def runko_default_handler():
    h = logging.StreamHandler()
    h.setLevel(logging.NOTSET)

    f = logging.Formatter('%(levelname)s(%(name)s): %(message)s')
    h.setFormatter(f)

    return h


if on_main_rank():
    runko_logger().setLevel(logging.INFO)

runko_logger().addHandler(runko_default_handler())
