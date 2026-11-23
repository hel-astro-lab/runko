# Copyright 2025 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Internal logging utilities which are configurable through standard logger module.

Logging in runko is done using named loggers `runko-N`,
where N is the rank where the logging comes from.

By default only main rank will have default logging level changed to INFO
while other ranks have logging level WARNING (python default).
The default handler used in runko is `logging.StreamHandler`
which logging level is set to `logging.NOTSET`.

Because the rank is not known before MPI is initialized,
the default configuration is applied lazily on the first call to `runko_logger`.
Determining the rank initializes the runko runtime (see `runko.runtime`).
"""


import logging
import sys
from .runtime import default_runtime


def on_main_rank() -> bool:
    """
    Checks if the caller is "main" rank.

    There is no other guarantees other that there is only one main rank.
    """

    return default_runtime().mpi_rank() == 0


def runko_logger_name(sublogger=None):
    name = f"runko-{default_runtime().mpi_rank()}"
    if sublogger:
        name += f".{sublogger}"
    return name


def runko_logger(sublogger=None):
    _configure_default_logging()
    return logging.getLogger(runko_logger_name(sublogger))


def runko_default_handler():
    h = logging.StreamHandler(stream=sys.stdout)
    h.setLevel(logging.NOTSET)

    f = logging.Formatter('%(levelname)s(%(name)s): %(message)s')
    h.setFormatter(f)

    return h


_default_logging_configured = False


def _configure_default_logging():
    global _default_logging_configured
    if _default_logging_configured:
        return
    _default_logging_configured = True

    logger = logging.getLogger(runko_logger_name())

    if on_main_rank():
        logger.setLevel(logging.INFO)

    logger.addHandler(runko_default_handler())
