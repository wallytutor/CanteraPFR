# CanteraPFR

Provides several variations of a plug-flow reactor model constructed over
Cantera framework of chemical kinetics.

# Installation

## Cygwin install (Windows)

**Note: this install comprises only C++ modules for now.**

Cygwin install is made easy by script `get_cygwin.bat`. This will download
locally everything you need to be available in order to compile the C++
executables. From you file navigator, just double click `get_cygwin.bat` and
wait for installation to proceed (this will download Cygwin and all required
packages, download Cantera and compile Cantera). The process lasts from 10 to
30 minutes typically, depending on your processor and internet connection.

Once the process has finished, double-click cygwin\\Cygwin.bat and just `cd`
(some NIX command line is required) to `CanteraPFR` directory. Once you have
reached your destination, run `make` and the executables must compile.

If you have no experience with Cygwin, take a small tutorial on how to move
across folders in Linux systems. *Under Cygwin, your drive C: is found under
/cygdrive/c and the same apply for other possible drives.*

Limitations:
- Currently Python modules are not compiled for Cygwin (TODO), so only XML
mechanism files will work with C++ executables. This will be fixed soon.

## Basic Python modules (Linux/Bash for Windows)

**Sorry, a conda build is not yet available, we are working on it!**

The installation assumes you are using Anaconda/Miniconda Python distribution.
For now the package is supported only on Nix operating systems because the
solver interface has not yet been adapted to easily compile under Windows (if
you are a *hardcore* programmer, you will know the steps to make it work under
native Windows, but for now I suggest you follow the step-by-step installation
provided in what follows).

If your platform is Windows 10, you need to enable Linux Subsystem for Windows
as described [here](https://docs.microsoft.com/en-us/windows/wsl/install-win10).
Select `Ubuntu` Linux distribution (or any other you feel comfortable with after
reading the next steps). Other distributions may work out-of-the-box but all
development and testing was carried with `Ubuntu`. Just keep in mind that here
we make use of `apt` for installing packages. The commands may change in other
distributions. Once that is done:

- Install [Miniconda](https://conda.io/miniconda.html) for your system.

Since current Python version is 3.7 and most probably this is the version your
Miniconda installer provided, you will need to create an environment with a
previous Python version (3.5). Here we will call this environment `py35`, but
you can call it whatever you want, just remember to replace `py35` by the name
you gave to your environment. This is necessary because one of the main packages
we use is outdated (we are working on a fix).

```
conda create -n py35 python=3.5
```

Before being able to use your Python 3.5, you need to activate the environment.
This is required to install the dependencies and every time you wish to use the
package.

```
conda activate py35
```

For some reason, Miniconda provides an old `pip` command. Upgrade it with:

```
pip install --upgrade pip
```

The wrapper to the solver we are going to use,
[scikits.odes](https://scikits-odes.readthedocs.io/en/latest/) requires a
specific version of [Sundials](https://computation.llnl.gov/projects/sundials)
library to work. Run the following command to check if you get version 2.7.0.
Otherwise, you will need to manually compile and install the library as given
in its documentation and make sure that libraries and header files are found
in your system's `PATH`/`LD_LIBRARY_PATH` variables. The documentation of
`scikits.odes` provides some guidelines about how to do this.

```
apt search libsundials-dev
```

If the version matches 2.7.0, you are good to go! Install Sundials with your
package manager.

```
sudo apt install libsundials-dev
```

With `conda` package manager install the main common dependencies:

```
conda install cython numpy scipy matplotlib sphinx pandas
```

Next install the core package [Cantera](https://cantera.org/):

```
conda install -c cantera cantera
```

Finally, using `pip` install the interface to the solvers:

```
pip install --user scikits.odes
```

Since the package is not complete, the `setup.py` script is not yet working.
You can place the whole package folder where your `PYTHONPATH` variable is able
to find it of append path from your user scripts (as the tests currently do).

## C++ executables

Compilation of C++ executables require main Cantera library to be available on
your system. Cantera is build with `scons`. If this tool is not available on
your system you can run:

```
pip install --user scons
```

or

```
conda install scons
```

If you already followed the steps for installing the Python modules, you already
have Sundials on your system. Otherwise you will need to install Sundials 2.7.0.
This is required because of the following bug:

**Compiling with Sundials > 30. Line 422 from `IDA_Solver.cpp`
`SUNLinSolFree((SUNLinearSolver) m_linsol);` leads to a segmentation fault.**

Now you are ready to run script `get_cantera.sh`. This script will download
and compile Cantera in your working directory (it will take a while). In order
to produce small executables, using shared libraries is recommended and you
need to append you `LD_LIBRARY_PATH` environmental variable with the full path
to `external/lib`. This can be done for your user with:

```
echo "export LD_LIBRARY_PATH=\"$(pwd)/external/lib:\$LD_LIBRARY_PATH\"" >> ~/.bashrc
```

or you can choose to manually edit the `Makefile` and use `-lcantera` instead
of `-lcantera_shared`.

# Documentation

In order to generate the documentation, go to directory `doc/` and run `make`.
This will provide a list of types you can compile the documentation. The default
and recommended format is `make html`. Documentation will be generated under
`doc/build/<type>`.

# Contribute

Please do not hesitate to contact me if you have suggestions or wish to request
new functionalities. I still have a lot of materials from my PhD to organize
and make available here. This will comprise, depending on my free time, a PFR
with imposed wall temperature, an adiabatic PFR (adapted from [here](https://github.com/Cantera/cantera-jupyter/blob/master/reactors/1D_pfr_surfchem.ipynb)), and an automatic
skeletal mechanism generation with DRG method.

Email: waltermateriais@gmail.com

# TODO

1. Follow [this](https://conda.io/docs/user-guide/tutorials/build-pkgs.html)
tutorial to get this package installed with `conda` one day.
1. Provide built-in DAE solver to avoid `scikits.odes`. Installation is limited
to Nix systems because it is difficult to provide a reliable way to install this
package that requires cython, BLAS/LAPACK, and Sundials to be available:
    - Notice that Cantera already has an interface to IDA. This could be
    cythonized and provide a quick solution to this. Compiling Sundials on
    Windows is simpler than doing a whole new solver.
