Documentation
=============

Documentation is generated with Sphinx and hosted on Read the Docs.

When modifying the documentation,
it is helpful to test the pages first on your own machine.
If CMake finds Sphinx at the configuration stage,
a build target ``docs`` will be available for building the documentation.

Sphinx and the required extensions can be installed with:

.. code:: shell

   pip install -r </path/to/runko>/docs/requirements.txt

Assuming the CMake build directory ``<build>`` has been configured with Sphinx,
the documentation can be built with:

.. code:: shell

   cmake --build <build> -t docs

The documentation should then be available at ``<build>/docs/docs/sphinx/index.html``.


Contributing to the Documentation
---------------------------------

Once the local documentation build is set up and working, adding new pages is easy:
drop new ``.rst`` files into the ``docs/`` folder and follow Sphinx's reStructuredText syntax.







