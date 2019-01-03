# CanteraPFR

Provides several variations of a plug-flow reactor (PFR) model constructed over
Cantera framework of chemical kinetics. These models range from isothermal to
predefined wall profiles variations of PFR.

# Installation

## Cygwin install (Windows)

**Note: this install comprises only C++ modules for now.**

If you do not have Cygwin:

Cygwin install is made easy by script `get_cygwin.bat`. This will download
locally everything you need to be available in order to compile the C++
executables. From you file navigator, just double click `get_cygwin.bat` and
wait for installation to proceed (this will download Cygwin and all required
packages, download Cantera and compile Cantera). The process lasts from 10 to
30 minutes typically, depending on your processor and internet connection.

Once the process has finished, double-click `cygwin\\Cygwin.bat` and just `cd`
(some NIX command line is required) to `CanteraPFR` directory. Once you have
reached your destination, run `make` and the library must compile out-of-the-box.
If you have no experience with Cygwin, take a small tutorial on how to move
across folders in Linux systems. *Under Cygwin, your drive C: is found under
/cygdrive/c and the same apply for other possible drives.*

If you have Cygwin:

Open `get_cygwin.bat` for edition and make sure all packages listed on line 22
of this file are available in your machine. Then execute your Cygwin terminal
and move to this folder. Manually run `get_cantera.sh` and wait for completion.
Once done, move to `CanteraPFR` and run `make` to generate library.

## Linux

Install `libopenblas-dev`, `libboost-dev` assuming you are using Ubuntu (or their
variants for your distro) and `scons` prior to installation. Next move to this
directory and run `get_cantera.sh` and wait for completion. Once done, move to
`CanteraPFR` and run `make` to generate the library.

## Note for Windows and Linux

You can choose to manually edit the `Makefile` and use `-lcantera` instead of
`-lcantera_shared`. This will produce larger executables and unless you bind
Sundials and OpenBLAS statically, provides no advantage except for the fact
that you will not need to append your `.bashrc` as automatically done by
`get_cantera.sh`. If you do not want to append your `.bashrc`, comment out
line 20 of `get_cantera.sh`.

## Create Python interface

Once the library was built, you can run the following command to get Python
modules installed:

```
python setup.py install
```

# Documentation

**Outdated/deprecated:** the documentation will continue to be as follows once
Cython modules are written/updated.

In order to generate the documentation, go to directory `doc/` and run `make`.
This will provide a list of types you can compile the documentation. The default
and recommended format is `make html`. Documentation will be generated under
`doc/build/<type>`.

# Contribute

Please do not hesitate to contact me if you have suggestions or wish to request
new functionalities. I still have a lot of materials from my PhD to organize
and make available here. This will comprise, depending on my free time, a PFR
with imposed wall temperature, an adiabatic PFR (adapted from
[here](https://github.com/Cantera/cantera-jupyter/blob/master/reactors/1D_pfr_surfchem.ipynb)),
and an automatic skeletal mechanism generation with DRG method.

Email: waltermateriais@gmail.com

# Limitations

1. Only XML mechanism files will work with C++ executables. This will be fixed.
1. Cantera shared library failed on Cygwin!
