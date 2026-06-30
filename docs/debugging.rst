Debugging
#########

Debuggers
=========

.. tabs::
   .. tab:: gdb
      .. code:: shell

         gdb [gdb-args...] --args <program> [args...]

      .. tip::
         gdb flag ``--tui`` can be useful.

      .. tip::
         gdb flag ``-ex <cmd>`` will run the given command at startup.
         Commonly it is convenient to use ``-ex run``.

   .. tab:: lldb

      .. code:: shell

         lldb [lldb-args...] -- <program> [args...]


Debuggin MPI programs
---------------------

One strategy for debugging MPI programs is to launch a new terminal for each MPI rank
and run debugger on each of those processes.

For example, if we are using OpenMPI (i.e. we are using ``mpirun``) and want to use xterm:

.. code:: shell

   mpirun -n <NP> xterm -hold -e <debugger command>

``-hold`` flag is useful to prevent the terminal process from exiting when closing the debugger.
``-e`` indicates the command to be executed in the new terminal.

See the documentation of your terminal for similar flags.
For example, with foot terminal:

.. code:: shell

   mpirun -n <NP> foot --hold <debugger command>


Valgrind
========

To use Valgrind with Python, you need to use
`Valgrind suppression file provided by the Python developers <https://github.com/python/cpython/blob/main/Misc/valgrind-python.supp>`_
or else there will be too many false positives due to Python's custom memory management.
Note that you have to uncomment lines for ``PyObject_Free`` and ``PyObject_Realloc`` in the suppression file
(`source <https://stackoverflow.com/questions/3982036/how-can-i-use-valgrind-with-python-c-extensions>`_).

.. code:: shell
   valgrind --tool=memcheck --suppressions=valgrind-python.supp python <script.py>


Valgrind + MPI
--------------

Similarly to using debuggers with MPI, we want to run Valgrind on each rank.
Its good idea to use log-file instead of screen to make sense of the output.

.. code:: shell

   mpirun -n <NP> valgrind --tool-memcheck --log-file=log.%p --suppressions=valgrind-python.supp python <script.py>
