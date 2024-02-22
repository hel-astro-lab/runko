Documentation
=============

Documentation is generated with ``sphinx`` and hosted at readthedocs. 


Local installation of Sphinx
----------------------------

In order to modify the documentation, it is helpful to test the pages first on your own machine. 

To do this, install sphinx with (on MacOS)

.. code-block:: bash

    brew install sphinx-docs


or in Ubuntu as


.. code-block:: bash

    apt-get install python3-sphinx

and finally the documentation theme


.. code-block:: bash
    
    pip3 install sphinx_rtd_theme --break-system-packages


The webpages can then be compiled (in the `build/` directory) as

.. code-block:: bash

   cmake ..
   make docs

After these steps, the documentations should be available in ``runko/build/docs/docs/sphinx/index.html``.


Contributing to the Documentation
---------------------------------

After you have the web page compilation setup up and running, creating new pages is easy: just add new .rst files into the ``docs/`` folder and follow the sphinx markdown syntax.







