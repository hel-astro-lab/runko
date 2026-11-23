# Copyright 2026 - 2026, Miro Palmu, Joonas Nättilä and the runko contributors
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Access to a process wide default runko runtime.

Importing runko does not initialize the runtime (nor MPI).
The runtime, and thus MPI, is initialized when the first
`runko.actions.RuntimeInstance` is created, either explicitly by the user
or implicitly by `default_runtime` when runko needs MPI information
(e.g. the rank for logging or tile partitioning).

All `RuntimeInstance` objects that are alive at the same time
share the same underlying runtime.
"""

import atexit
import runko_cpp_bindings.actions as actions

_default_runtime = None


def _release_default_runtime():
    global _default_runtime
    _default_runtime = None


def default_runtime() -> actions.RuntimeInstance:
    """
    Returns the process wide default runtime instance,
    creating it (and initializing MPI) on the first call.

    The default instance is released at interpreter exit.
    """

    global _default_runtime
    if _default_runtime is None:
        _default_runtime = actions.RuntimeInstance()
        atexit.register(_release_default_runtime)
    return _default_runtime
