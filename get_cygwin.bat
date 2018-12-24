@echo off

:: ----------------------------------------------------------------------------
:: Download Cygwin installer
:: ----------------------------------------------------------------------------

PowerShell.exe -ExecutionPolicy Bypass -Command "(New-Object Net.WebClient).DownloadFile('https://cygwin.com/setup-x86_64.exe', 'setup-x86_64.exe')"

:: ----------------------------------------------------------------------------
:: Install Cygwin with all dependencies
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
    --packages make,wget,gcc-core,gcc-g++,libopenblas,libsundials-devel,libboost-devel,python3,python3-devel,python3-cython,python3-numpy,scons,gnuplot

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
