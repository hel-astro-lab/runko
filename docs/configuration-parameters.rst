.. _conf-params:

Configuration Parameters
########################

.. role:: python(code)
   :language: python

Simulations require a lot of configuration parameters to be set (e.g. mesh size, what algorithms to use, etc.).
In Runko, these are set using :python:`runko.Configuration` Python class.
It can be thought as a Python dictionary but the values are accessed through attribute syntax
instead of the dictionary syntax (:python:`config.key = value` vs. :python:`some_dict['key'] = value`).
If you want to access a value using a string as the key, you can use
`getattr <https://docs.python.org/3.11/library/functions.html#getattr>`_ function.

Different parts of runko take :python:`runko.Configuration` objects as arguments.
Runko will raise an error if a given configuration object does not contain the required parameters.
Parameters can also be optional. Missing a optional parameter is not a error,
but attempting to use functionality that depends on that parameter will result in an error.
This can be annoying if it happens at the end of a large simulation,
so make sure to always test your runs on a small scale first.


Constructing :python:`runko.Configuration` object
=================================================

:python:`runko.Configuration` can be constructed completely inline in Python,
which is some times useful for 100% self-contained simulation files:

.. code:: python

   import runko
   config = runko.Configuration(None)
   config.foo = 42
   config.bar = "abc"
   # ...


More commonly, the configuration object is constructed from a configuration file
by passing the file path to the constructor.
The configuration file syntax loosely follows what Python :python:`configparser` supports by default.
In short, there is one key–value pair per line, and the key and value are separated by either
``:`` or ``=``. Comments start with either ``#`` or ``;``.
Key-value pairs are required to be organized into sections which are marked using ``[section-name]``.
However, the sections are completely optional and have no semantic meaning.
Starting from Python >= 3.13, sections are made completely optional.
A common convention is to have following sections: io, simulation, grid, problem, particles and algorithms.

Here is an example of a configuration file:

.. code:: ini

   [io]
   # This is a comment.
   ; And this too.
   foo: 42
   bar = "abc"

   [simulation]
   # ...


And here is an example of how to read it:

.. code:: python

   import runko
   config = runko.Configuration("path/to/config")

   # Sometimes it is useful to set other parameters dynamically
   # based on other parameters.

   config.foobar = config.foo * config.bar


Supported values
================

Internally, Runko interprets the values from the configuration file as Python literals,
and therefore all Python literal data types are supported.
However, different parts of Runko expect parameters to be of specific types
and will raise an error if the expected type does not match the actual type.
For example, if Runko expects an integer parameter, then the value ``2.0`` will results in an error.


Inspecting configuration from output directory
==============================================

For each simulation, runko will copy the configuration file to the output directory if it exists.
In any case, runko will write a `pickled <https://docs.python.org/3/library/pickle.html>`_
version of the configuration object to the output directory as ``<output-dir>/config.pkl``.

Contents of the pickled configuration object can be inspected using the ``runko`` tool:

.. code:: shell

   runko inspect conf <path/to/output/dir>


List of common parameters
=========================

As mentioned above, there are different types of parameters in Runko.
Some parameters are always required, while others are optional.
Additionally, some parameters are not required by any internal Runko machinery
but are still useful
(e.g. parameters used during the simulation, for post‑processing, or for identifying different simulations).
If "used in"-field is empty, the parameter is not actually used by runko,
but its use in the project files has become convention.

.. csv-table:: IO parameters
   :header: "name", "data type", "used in", "description"
   :file: params/io.csv

.. csv-table:: Simulation parameters
   :header: "name", "data type", "used in", "description"
   :file: params/sim.csv

.. csv-table:: Miscellaneous parameters
   :header: "name", "data type", "used in", "description"
   :file: params/misc.csv

.. csv-table:: Algorithm parameters
   :header: "name", "data type", "used in", "description"
   :file: params/algs.csv
