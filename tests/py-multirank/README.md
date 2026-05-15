This directory contains tests which require multiple MPI ranks.

`mpi_unittest.py` contains utility functions for multirank assertions.

# Running tests

Inside configured CMake build directory run: `ctest`

# Adding tests
## Python test

In `CMakeLists.txt` add `add_python_mpi_test` call.
It takes two arguments:

- name of the test
  - corresponding test file needs to be in form: `test_<test_name>.py`
- how many tasks should be used to run the test


