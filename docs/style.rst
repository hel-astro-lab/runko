Coding style
============

Runko tries to conform to common C++ and Python3 styles. 
The motivation is to provide 

i. a uniform code base that is familiar to everybody, 
ii. enforce a coding style that is less error prone, and 
iii. encourage more succint notation since less code means less bugs.

These styles can be enforced and checked with various tools, that are run manually from time to time.


General guidlines
-----------------

* use "snake_case" for variable and function names, i.e., ``ex``, ``gam``, ``get_container()``
* objects and classes are named with a Capital letter and using the CamelCase form, i.e., ``Tile``, ``ParticleContainer``
* use short but descriptive naming that resembles the underlying mathematical formulas 
* stick to <100 character long lines. 
* surround ``+`` and ``-`` operators with a space for easier visibility, i.e., ``2*x + y``.
* do **not** use tabs but normal spaces for indentation.
* Python indentation level is ``4`` spaces and C++ is ``2``; this enforces more compact code.



Black
-----

`black` is an uncompromised Python code formatter that is also PEP 8 compliant.
In general, we try to conform to the output of ``black`` as much as possible, but we do deviate sometimes, in the name of more succint expressions. 

Format a python script by running

.. code-block:: bash:

    black some.py

.. note::
    Black overwrites the file. It is typically easiest to apply it by first developing your code as normal, then commiting to the git. After that, finish by running black on your python files. You can then check the changes with ``git diff`` and see that they make sense.




clang-tools
-----------

Clang-tools are a set of code formatting, analysis, and testing tools based on the LLVM compiler set that aim to make C++ code more readable and less bugprone.

In order to use the tools, first we need to create a ``compile-commands.json`` file that provides a recipe for the LLVM compiler how to compile each file.
For details, see [CD1]_, [CD2]_.

To generate the ``compile-commands.json`` with CMake, run it in the ``runko/build/`` directory with

.. code-block:: bash:

    cmake -DCMAKE_EXPORT_COMPILE_COMMANDS=1 ..

And prune it manually by removing unwanted files (like those of the included libraries, ``corgi/*``, ``corgi/mpi4cpp/`` etc.).



Clang-Format
^^^^^^^^^^^^


The `C++` files are regularly checked and formatted with `Clang-Format <https://clang.llvm.org/docs/ClangFormat.html>`_.
Runko ships with a ``clang-format`` configuration file at the top level of the repository (filename ``.clang-format``).

Currently formatting is **not** applied automatically, but manually using ``clang-format`` for newly developed files is highly encouraged.
To check if a file needs formatting,

.. code-block:: bash:

    clang-format -style=file --dry-run some.c++

The output will show what needs to be changed, if any. To automatically apply the changes,


.. code-block:: bash:

    clang-format -style=file -i some.c++


Clang-Tidy
^^^^^^^^^^

`Clang-Tidy <https://clang.llvm.org/extra/clang-tidy/>`_ is a similar but more heavy-weight tool that performs a deeper static code analyses and checks for bugprone code etc.

We typically test the code with checks from:

* ``modernize-*``
* ``performance-*``
* ``bugprone-*``
* ``readability-*``
* ``mpi-*``
* ``cppcoreguidlines-*``

in roughly the given order.

A helper python script exists in `tools/` to perform the analysis for all files. Run the script (still in the ``runko/build/`` directory) with


.. code-block:: bash:

    python3 ../tools/run-clang-tidy.py -header-filter='.*' -checks='-*,modernize-*'
    python3 ../tools/run-clang-tidy.py -header-filter='.*' -checks='-*,performance-*'
    python3 ../tools/run-clang-tidy.py -header-filter='.*' -checks='-*,readability-*'
    python3 ../tools/run-clang-tidy.py -header-filter='.*' -checks='-*,bugprone-*'
    python3 ../tools/run-clang-tidy.py -header-filter='.*' -checks='-*,cppcoreguidlines-*'

Use common sense what fixes to actually apply.



Include-what-you-use
--------------------

To run include what you use, install (``brew install include-what-you-use`` on macOS), then run:

.. code-block:: bash:

    cmake -S . -B build-iwyu -DCMAKE_CXX_INCLUDE_WHAT_YOU_USE=$(which include-what-you-use)

    cmake --build build

.. warning::
    TODO: not implemented yet.








References
----------

.. [CD1] `<https://sarcasm.github.io/notes/dev/compilation-database.html>`_

.. [CD2] `<https://eli.thegreenplace.net/2014/05/21/compilation-databases-for-clang-based-tools>`_




