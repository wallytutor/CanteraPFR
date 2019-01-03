###############################################################################
##  MAKEFILE
###############################################################################

include Makefile.in

SOURCES   := src
# SUNDIALS  := -lsundials_ida -lsundials_nvecserial
CANTERA   := -L$(EXTERNAL)/lib -lcantera_shared
INCLUDES  := -I$(EXTERNAL)/include -Iinclude
LIBRARIES := $(CANTERA) -lopenblas -pthread
CPPFLAGS  := $(OPTIONS) -fPIC $(INCLUDES)
CXX       := g++
AR        := ar
ARFLAGS   := sq

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

test: $(BINARIES)

xAdiabaticPFR: $(OBJECTS)
	$(CXX) $(CPPFLAGS) test/main/xAdiabaticPFR.cpp src/AdiabaticPFR.o \
	$(LIBRARIES) -o bin/xAdiabaticPFR

xHeatWallPFR: $(OBJECTS)
	$(CXX) $(CPPFLAGS) test/main/xHeatWallPFR.cpp src/HeatWallPFR.o \
	$(LIBRARIES) -o bin/xHeatWallPFR

xIsothermalPFR: $(OBJECTS)
	$(CXX) $(CPPFLAGS) test/main/xIsothermalPFR.cpp src/IsothermalPFR.o \
	$(LIBRARIES) -o bin/xIsothermalPFR

clean:
	$(RM) $(OBJECTS)

dist-clean: clean
	$(RM) bin/* lib/*
