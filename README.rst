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

Derived Settings
----------------
The following settings are derived from other settings and are not needed in
the configuration file:

.. code-block:: ini

    [&system]
    nat = [ATOMIC_POSITIONS]
    ntyp = [ATOMIC_SPECIES]

    [sigma]
    band_index_min = [pp_in].vxc_diag_nmin
    band_index_max = [pp_in].vxc_diag_nmax

    [absorption]
    number_val_bands_coarse = [kernel].number_val_bands
    number_cond_bands_coarse = [kernel].number_cond_bands
