mpi4py might have to be installed manually.
The default installation, obtained automatically or
with ``python -m pip install mpi4py``,
can link against a different MPI implementation than the one runko uses.
This can lead to problems, so mpi4py can instead be installed manually:

.. code:: shell

   MPI4PY_BUILD_MPICC=<mpi-aware C compiler> python -m pip install --no-cache-dir --no-binary=mpi4py mpi4py
