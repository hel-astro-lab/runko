Developer workflow
##################

Building
========

The basic :ref:`user installation of runko <building>`
is not sufficient when developing new features or when running tests.
There are two other options which give developer more control of the build process.
The first and recommended way is editable install with pip.
The second is manual CMake usage which allows to run tests
but does not allow to installing the runko python package.
Editable installs allows to do both.

For both of the methods runko repository has to be manually obtained:

.. code:: shell

   git clone https://github.com/hel-astro-lab/runko --recurse-submodules

Editable install
----------------

Editable installs allows to modify the python source files of runko
such that the modifications take effect without having to install runko again.
However, if the C++ source files are modified, then runko has to be installed again,
but the difference is that the CMake build directory can be reused
and only the modified parts have to compiled (assuming ``--no-build-isolation``).
The persistent CMake build directory also allows to run tests.

The editable builds with ``--no-build-isolation`` requires manual installation
of the build time dependencies of the python package:

.. code:: shell

   pip install pybind11 scikit-build-core

Runko ships per-machine configure presets in ``CMakePresets.json``
(``unix-cpu``, ``macos-cpu``, ``unix-gpu``, ``lumi-cpu``, ``lumi-gpu``,
``hile-cpu``, ``hile-gpu``); pick the one matching your machine and
backend, and forward it to CMake via ``cmake.args``. With hip backend on
a generic Linux box as an example:

.. code:: shell

   pip install --no-build-isolation -v -e <path/to/runko/repo> \
     --config-settings=cmake.args=--preset=unix-gpu


``--no-build-isolation`` means that the build directory is persistent between the installs,
and only the modified parts of the C++ has to compiled again on the next install.
Verbose flag ``-v`` is useful to see full outputs from CMake.

Preset cache variables can be overridden per invocation by appending
``--config-settings=cmake.define.<DEFINE>=<VALUE>`` or
``--config-settings=cmake.build-type=<BUILD-TYPE>`` — the override always
wins, so the same preset can be flipped to a Debug build:

.. code:: shell

   pip install --no-build-isolation -v -e <path/to/runko/repo> \
     --config-settings=cmake.args=--preset=unix-gpu \
     --config-settings=cmake.build-type=Debug

The scikit-build-core managed build directory is
``<path/to/runko/repo>/build/{state}-{build_type}-{wheel_tag}`` and is
reused between editable installs.

For more configuration options, see `scikit-build-core docs <https://scikit-build-core.readthedocs.io/en/latest/>`_.


Manual CMake usage
------------------

Manual CMake usage will enable full functionality of CMake,
but it does not enable installing the runko python package.
Use a preset to configure and build — the preset name is also the build
directory name (each preset writes to ``<repo>/<preset-name>/`` via the
shared ``_common`` hidden base):

.. code:: shell

   cmake --preset=unix-gpu -S <path/to/runko/repo>
   cmake --build <path/to/runko/repo>/unix-gpu

Override the build type or any cache variable with ``-D`` on either
command, e.g. ``-DCMAKE_BUILD_TYPE=Debug``.


Testing
=======

Running tests
-------------

Runko uses `CMake's ctest <https://cmake.org/cmake/help/latest/manual/ctest.1.html>`_ to run tests.
Once a runko build directory has been configured, either with editable install or manually with CMake,
runko tests can be run with ``ctest`` inside the build directory.

See the ctest documentation for all the options, but here is some of the most useful ones:

================== ==========================================
Flag               Description
================== ==========================================
``-j N``           Run N tests in parallel
------------------ ------------------------------------------
``-R <regex>``     Run tests matching regular expression.
------------------ ------------------------------------------
``--verbose``      Enable verbose output from tests.
------------------ ------------------------------------------
``--rerun-failed`` Run only the tests that failed previously.
================== ==========================================
