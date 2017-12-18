.. image:: https://travis-ci.org/nelsyeung/bgwtools.svg?branch=master
    :target: https://travis-ci.org/nelsyeung/bgwtools
.. image:: https://coveralls.io/repos/github/nelsyeung/bgwtools/badge.svg
    :target: https://coveralls.io/github/nelsyeung/bgwtools

BerkeleyGW Tools
================
Generate BerkeleyGW and QuantumESPRESSO simulation files from a single input

Installation
------------
.. code-block:: bash

    git clone https://github.com/nelsyeung/bgwtools.git
    pip install -e .

Usage
-----
Create an bgwtools.ini input file. Check examples/bgwtools.ini.

.. code-block:: bash

    bgwtools bgwtools.ini

This will create a bgw and a qe folder with all the required input files.
