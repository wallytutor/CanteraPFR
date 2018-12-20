# CanteraPFR

Provides several variations of a plug-flow reactor model constructed over
Cantera framework of chemical kinetics.

# Installation

**Sorry, a conda build is not yet available, we are working on it!**

The installation assumes you are using Anaconda/Miniconda Python distribution.
For now the package is supported only on Nix operating systems because the
solver interface has not yet been adapted to easily compile under Windows (if
you are a *hardcore* programmer, you will know the steps to make it work under
native Windows, but for now I suggest you follow the step-by-step intallation
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
conda install cython numpy scipy matplotlib sphinx
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

# TODO

1. Provide a full C++ implementation of basic calculators.
1. Follow [this](https://conda.io/docs/user-guide/tutorials/build-pkgs.html)
tutorial to get this package installed with `conda` one day.
1. Provide built-in DAE solver to avoid `scikits.odes`. Installation is limited
to Nix systems because it is difficult to provide a reliable way to install this
package that requires cython, BLAS/LAPACK, and Sundials to be available.
