mpi4py might have to be installed manually.
The default installation obtained automatically or
with ``python -m pip install mpi4py``
can link to different mpi implementation compared what runko is using.
This can lead to problems, so mpi4py can be installed manually:

.. code:: shell

   MPI4PY_BUILD_MPICC=<mpi-aware C compiler> python -m pip install --no-cache-dir --no-binary=mpi4py mpi4py
