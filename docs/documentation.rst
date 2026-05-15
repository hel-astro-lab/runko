Documentation
=============

Documentation is generated with Sphinx and hosted at readthedocs.

In order to modify the documentation,
it is helpful to test the pages first on your own machine.
If CMake finds Sphinx at configuration stage
then there will be a build target ``docs`` which will build the documentation.

Spinx and required extensions can be installed with:

.. code:: shell

   pip install -r </path/to/runko>/docs/requirements.txt

Assuming CMake build directory ``<build>`` has been configured with Sphinx,
documentation can be build with:

.. code:: shell

   cmake --build <build> -t docs

The documentations should be available at: ``<build>/docs/docs/sphinx/index.html``


Contributing to the Documentation
---------------------------------

After you have the web page compilation setup up and running, creating new pages is easy:
just add new .rst files into the ``docs/`` folder and follow the sphinx markdown syntax.







