###############################################################################
##  MAKEFILE
###############################################################################

PYTHON       := python
OPTIONS      := -g -O3 -std=c++11
BLAS_LAPACK  := -lopenblas
LIB_SUNDIALS := -Lexternal/Nix/lib -lsundials_ida -lsundials_nvecserial
INC_SUNDIALS := -Iexternal/Nix/include
LIB_CANTERA  := -Lexternal/Nix/lib -lcantera_shared
INC_CANTERA  := -Iexternal/Nix/include
LIB_SELF     := -LCanteraPFR/lib
INC_SELF     := -ICanteraPFR/include
LIBRARIES    := $(LIB_CANTERA) $(LIB_SUNDIALS) $(LIB_SELF) $(BLAS_LAPACK)
INCLUDES     := $(INC_CANTERA) $(INC_SUNDIALS) $(INC_SELF)


LIBRARIES += -pthread
CPPFLAGS  := $(OPTIONS) -fPIC $(INCLUDES)
SOURCES   := CanteraPFR/src/cpp

CXX       := g++
AR        := ar
ARFLAGS   := sq
RM        := rm -rf

###############################################################################
##  CONVERT OBJECTS
###############################################################################

SOBJECT  := $(shell find $(SOURCES) -iname "*.cpp")
OBJECTS  := $(patsubst $(SOURCES)/%.cpp,$(SOURCES)/%.o,$(SOBJECT))

###############################################################################
##  COMPILE MODELS
###############################################################################

.PHONY := all
all: cython docs shared_library

cython: library
	$(PYTHON) setup.py install

docs: cython
	cd docs && make html

library: $(OBJECTS)
	$(AR) $(ARFLAGS) CanteraPFR/lib/libCanteraPFR.a $(OBJECTS)

shared_library: $(OBJECTS)
	$(CXX) -shared -fPIC $(INCLUDES) CanteraPFR/cPFR.c $(OBJECTS) $(LIBRARIES) \
	-o CanteraPFR/lib/libCanteraPFR_shared.so

clean:
	$(RM) $(OBJECTS) build *.egg-info

dist-clean: clean
	$(RM) CanteraPFR/lib/*
	cd docs && make clean

clean-all: dist-clean
	$(RM) external cygwin setup-x86_64.exe

###############################################################################
##  EOF
###############################################################################
