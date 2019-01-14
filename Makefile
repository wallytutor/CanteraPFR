###############################################################################
##  MAKEFILE
###############################################################################

ifeq (,$(wildcard Makefile.in))
include Makefile.in.default
else
include Makefile.in
endif

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
