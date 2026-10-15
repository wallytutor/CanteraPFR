@echo off

:: ----------------------------------------------------------------------------
:: Download Cygwin installer
:: ----------------------------------------------------------------------------

PowerShell.exe^
    -ExecutionPolicy Bypass^
    -Command "(New-Object Net.WebClient).DownloadFile('https://cygwin.com/setup-x86_64.exe', 'setup-x86_64.exe')"

:: ----------------------------------------------------------------------------
:: Install Cygwin with all dependencies
:: https://www.scivision.co/matplotlib-in-cygwin-64-bit/ -> matplotlib
:: ----------------------------------------------------------------------------

setup-x86_64.exe^
    --no-admin^
    --no-desktop^
    --no-shortcuts^
    --no-startmenu^
    --only-site^
    --site http://cygwin.mirror.constant.com^
    --quiet-mode^
    --root cygwin^
    --packages vim,git,make,wget,gcc-core,gcc-g++,gcc-fortran^
    --packages libboost-devel,libopenblas,libsundials-devel^
    --packages python3,python3-devel,python3-pip,python3-cython^
    --packages python3-setuptools,python3-numpy,python3-sphinx,gnuplot^
    --packages pkg-config,ghostscript,libfreetype-devel,libpng-devel^
    --packages libgtk2.0-devel,openbox,python3-pyqt5,cmake,cmake-gui

:: ----------------------------------------------------------------------------
:: This script's path.
:: ----------------------------------------------------------------------------

set DRIVE=%~d0
set CURRENTDIR=%~p0
set CURRENTDIR=%CURRENTDIR:\=/%

:: ----------------------------------------------------------------------------
:: Cygwin path of the current directory.
:: ----------------------------------------------------------------------------

set THISRUNDIR="/cygdrive/%DRIVE:~0,-1%%CURRENTDIR%"

:: ----------------------------------------------------------------------------
:: Retrieve and compile Cantera.
:: ----------------------------------------------------------------------------

cygwin\bin\bash --login -c  "cd %THISRUNDIR%; ./get_cantera.sh;"

pause
