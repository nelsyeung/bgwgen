.. image:: https://travis-ci.org/nelsyeung/bgwgen.svg?branch=master
    :target: https://travis-ci.org/nelsyeung/bgwgen
.. image:: https://coveralls.io/repos/github/nelsyeung/bgwgen/badge.svg
    :target: https://coveralls.io/github/nelsyeung/bgwgen

BerkeleyGW Generator
====================
Generate BerkeleyGW and QuantumESPRESSO simulation files from a single input

Installation
------------
.. code-block:: bash

    git clone https://github.com/nelsyeung/bgwgen.git
    pip install -e .

Usage
-----
Create an bgwgen.ini input file. Check examples/bgwgen.ini.

.. code-block:: bash

    bgwgen bgwgen.ini

This will create a bgw and a qe folder with all the required input files.
