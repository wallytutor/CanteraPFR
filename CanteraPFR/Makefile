###############################################################################
##  MAKEFILE
###############################################################################

include Makefile.in

SOURCES   := src/cpp
SUNDIALS  := -lsundials_ida -lsundials_nvecserial
CANTERA   := -L$(EXTERNAL)/lib $(SUNDIALS) -lcantera_shared
INCLUDES  := -I$(EXTERNAL)/include -Iinclude
LIBRARIES := $(CANTERA) -lopenblas -pthread
CPPFLAGS  := $(OPTIONS) -fPIC $(INCLUDES)
CXX       := g++
AR        := ar
ARFLAGS   := sq
RM        := rm -rf

###############################################################################
##  CONVERT OBJECTS
###############################################################################

SOBJECT  := $(shell find $(SOURCES) -iname "*.cpp")
OBJECTS  := $(patsubst $(SOURCES)/%.cpp,$(SOURCES)/%.o,$(SOBJECT))
BINARIES := xAdiabaticPFR xHeatWallPFR xIsothermalPFR

###############################################################################
##  COMPILE MODELS
###############################################################################

.PHONY := all
all: $(OBJECTS)
	$(AR) $(ARFLAGS) lib/libCanteraPFR.a $(OBJECTS)
	python setup.py install && cd doc; make html

test: $(BINARIES)

xAdiabaticPFR: $(OBJECTS)
	$(CXX) $(CPPFLAGS) test/main/xAdiabaticPFR.cpp src/cpp/AdiabaticPFR.o \
	$(LIBRARIES) -o bin/xAdiabaticPFR

xHeatWallPFR: $(OBJECTS)
	$(CXX) $(CPPFLAGS) test/main/xHeatWallPFR.cpp src/cpp/HeatWallPFR.o \
	$(LIBRARIES) -o bin/xHeatWallPFR

xIsothermalPFR: $(OBJECTS)
	$(CXX) $(CPPFLAGS) test/main/xIsothermalPFR.cpp src/cpp/IsothermalPFR.o \
	$(LIBRARIES) -o bin/xIsothermalPFR

clean:
	$(RM) $(OBJECTS) build

dist-clean: clean
	$(RM) bin/* lib/*
