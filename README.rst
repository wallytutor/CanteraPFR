CanteraPFR
==========

Provides several variations of a plug-flow reactor (PFR) model constructed over
Cantera framework of chemical kinetics. These models range from isothermal to
predefined wall profiles variations of PFR.

Installation
============

For detailed instructions, read `INSTALL`_.

.. INSTALL: https://github.com/waltermateriais/CanteraPFR/blob/master/INSTALL.rst

Linux
-----

Tested on Ubuntu 18.04.

1. Install `libopenblas-dev`, `libboost-dev` (or equivalent) with your package
manager.
1. Install `scons` with `pip`, `conda` or your package manager.
1. Run `./get_cantera.sh` and wait for completion.
1. Run `make` and wait for completion.

Cygwin (Windows)
----------------

Tested on Windows 10.

1. Double-click `get_cygwin.bat` and wait for completion.
1. Double-click `cygwin\\Cygwin.bat` to open a terminal.
1. In terminal, move to this directory with `cd`.
1. Run `make` and wait for completion.

Documentation
=============

**Only after building the project, otherwise Sphinx will not be able to import!**

In order to generate the documentation, go to directory `CanteraPFR/doc/` and
run `make`. This will provide a list of types you can compile the documentation.
The default format is `make html` and this is already compiled with library is
built. Documentation will be generated under `CanteraPFR/doc/build/<type>`.

Contribute
==========

Please do not hesitate to contact me if you have suggestions or wish to request
new functionalities. I still have a lot of materials from my PhD to organize
and make available here. This will comprise, depending on my free time, a PFR
with imposed wall temperature, an adiabatic PFR (adapted from
[here](https://github.com/Cantera/cantera-jupyter/blob/master/reactors/1D_pfr_surfchem.ipynb)),
and an automatic skeletal mechanism generation with DRG method.

Email: waltermateriais@gmail.com

Limitations
===========

1. Cantera shared library failed on Cygwin!
1. For some weird reason Sundials 3.1.0 does not work with Cython modules (although
binaries compiled with it run just fine).
