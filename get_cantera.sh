#!/usr/bin/env bash

export EXTERNAL="$(pwd)/external"

git clone --recursive https://github.com/Cantera/cantera.git
cd cantera

git checkout tags/v2.4.0
git submodule update

scons build \
    prefix=$EXTERNAL \
    python_package="none" \
    system_sundials="y" \
    cxx_flags="-g -Wextra -O3 -std=c++11"

scons test
scons install
