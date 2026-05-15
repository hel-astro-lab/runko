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

Now to runko (with hip backend as an example) can be installed in editable mode with:

.. code:: shell

   pip install --no-build-isolation -v -e <path/to/runko/repo> \
   --config-settings=cmake.define.tyvi_BACKEND=hip \
   --config-settings=cmake.define.CMAKE_CXX_COMPILER=hipcc \
   --config-settings=cmake.build-type=Debug


``--no-build-isolation`` means that the build directory is persistent between the installs,
and only the modified parts of the C++ has to compiled again on the next install.
Verbose flag ``-v`` is useful to see full outputs from CMake.

As seen from above, we can give custom defines to CMake with ``--config-settings=cmake.define.<DEFINE>=<VALUE>``,
and set the CMake build type with ``--config-settings=cmake.build-type=<BUILD-TYPE>``.
The created build directory will be of form ``<path/to/runko/repo>/build/{state}-{build_type}-{wheel_tag}``.

For more configuration options, see `scikit-build-core docs <https://scikit-build-core.readthedocs.io/en/latest/>`_.

.. note::

   TODO: check editable installs work on macos??


Manual CMake usage
------------------

Manual CMake usage will enable full functionality of CMake,
but it does not enable installing the runko python package.
Documenting CMake usage is out of scope of this documentation.
Here is a simple configuration for hip backend:

.. code:: shell

   cmake -B debug-build -D tyvi_BACKEND=hip -D CMAKE_CXX_COMPILER=hipcc -D CMAKE_BUILD_TYPE=Debug <path/to/runko/repo>


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
