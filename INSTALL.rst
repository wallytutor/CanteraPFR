Install
=======

The following steps assume you do not dispose of a Cantera (at least 2.4.0 install in
your machine. **Notice that the following will only work if you compiled Cantera
locally with Sundials 2.7.0**. If you already have Cantera, just create a file named
`Makefile.in` in this same directory (same level as `Makefile`) and write the
following contents to this file:

```
EXTERNAL := full path to your Cantera and Sundials-2.7.0 roots (must be the same).
OPTIONS  := -g -O3 -std=c++11
```

If Sundials is elsewhere, you will need to edit the `Makefile` itself, and if you
intend to do this I assume you know how to do that.  If compiling from Cygwin, also
add `-U__STRICT_ANSI__` to `OPTIONS` line.

Linux
-----

Install `libopenblas-dev`, `libboost-dev` assuming you are using Ubuntu (or their
variants for your distro) and `scons` prior to installation. Next move to this
directory and run `get_cantera.sh` and wait for completion. Once done, run
`make` to generate the library.

Cygwin (Windows)
----------------

**Note: during development Cygwin is not supported. If you really need this
interface for now, contact me directly so that we can manage to install it.**

If you do not have Cygwin:

Cygwin install is made easy by script `get_cygwin.bat`. This will download
locally everything you need to be available in order to compile the C++
executables. From you file navigator, just double click `get_cygwin.bat` and
wait for installation to proceed (this will download Cygwin and all required
packages, download Cantera and compile Cantera). The process lasts from 10 to
30 minutes typically, depending on your processor and internet connection.

Once the process has finished, double-click `cygwin\\Cygwin.bat` and wait for
completion. Run `make` and the library must compile out-of-the-box.
If you have no experience with Cygwin, take a small tutorial on how to move
across folders in Linux systems. *Under Cygwin, your drive C: is found under
/cygdrive/c and the same apply for other possible drives.*

If you have Cygwin:

Open `get_cygwin.bat` for edition and make sure all packages listed on line 22
of this file are available in your machine. Then execute your Cygwin terminal
and move to this folder. Manually run `get_cantera.sh` and wait for completion.
Once done, run `make` to generate library.

Note for Cygwin and Linux
-------------------------

You can choose to manually edit the `Makefile` and use `-lcantera` instead of
`-lcantera_shared`. This will produce larger executables and unless you bind
Sundials and OpenBLAS statically, provides no advantage except for the fact
that you will not need to append your `.bashrc` as automatically done by
`get_cantera.sh`. If you do not want to append your `.bashrc`, comment out
line 20 of `get_cantera.sh`.
