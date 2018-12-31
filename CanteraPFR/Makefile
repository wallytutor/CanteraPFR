###############################################################################
##  MAKEFILE
###############################################################################

include Makefile.in

BINARIES  := xAdiabaticPFR \
						 xHeatWallPFR \
					   xIsothermalPFR

SOURCES   := src
SUNDIALS  := -lopenblas -lsundials_ida -lsundials_nvecserial
CANTERA   := -L$(EXTERNAL)/lib -lcantera_shared
INCLUDES  := -I$(EXTERNAL)/include -Iinclude
LIBRARIES := $(CANTERA) $(SUNDIALS) -pthread
CPPFLAGS  := $(OPTIONS) $(INCLUDES)
CXX       := g++

###############################################################################
##  CONVERT OBJECTS
###############################################################################

SOBJECT := $(shell find $(SOURCES) -iname "*.cpp")
OBJECTS := $(patsubst $(SOURCES)/%.cpp,$(SOURCES)/%.o,$(SOBJECT))

###############################################################################
##  COMPILE MODELS
###############################################################################

.PHONY := all
all: $(BINARIES)

xAdiabaticPFR: $(OBJECTS) src/AdiabaticPFR/*
	$(CXX) $(CPPFLAGS) src/AdiabaticPFR/*.o $(LIBRARIES) -o xAdiabaticPFR

xHeatWallPFR: $(OBJECTS) src/HeatWallPFR/*
	$(CXX) $(CPPFLAGS) src/HeatWallPFR/*.o $(LIBRARIES) -o xHeatWallPFR

xIsothermalPFR: $(OBJECTS) src/IsothermalPFR/*
	$(CXX) $(CPPFLAGS) src/IsothermalPFR/*.o $(LIBRARIES) -o xIsothermalPFR

clean:
	$(RM) $(OBJECTS)

dist-clean: clean
	$(RM) $(BINARIES)
