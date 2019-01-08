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
all: $(OBJECTS)
	$(AR) $(ARFLAGS) CanteraPFR/lib/libCanteraPFR.a $(OBJECTS)
	$(PYTHON) setup.py install && cd doc; make html

clean:
	$(RM) $(OBJECTS) build CanteraPFR.egg-info dist

dist-clean: clean
	$(RM) CanteraPFR/lib/*
	cd doc && make clean

clean-all: dist-clean
	$(RM) external cygwin setup-x86_64.exe

###############################################################################
##  EOF
###############################################################################
