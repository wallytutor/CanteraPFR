###############################################################################
##  MAKEFILE
###############################################################################

include Makefile.in

SOURCES   := CanteraPFR/src/cpp
SUNDIALS  := -lsundials_ida -lsundials_nvecserial
CANTERA   := -L$(EXTERNAL)/lib $(SUNDIALS) -lcantera_shared
INCLUDES  := -I$(EXTERNAL)/include -ICanteraPFR/include
LIBRARIES := $(CANTERA) -lopenblas -pthread
CPPFLAGS  := $(OPTIONS) -fPIC $(INCLUDES)
CXX       := g++
AR        := ar
ARFLAGS   := sq
RM        := rm -rf

TEST      := CanteraPFR/test/main
BINS      := CanteraPFR/bin

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
	$(AR) $(ARFLAGS) CanteraPFR/lib/libCanteraPFR.a $(OBJECTS)
	python setup.py install && cd CanteraPFR/doc; make html

test: $(BINARIES)

xAdiabaticPFR: $(OBJECTS)
	$(CXX) $(CPPFLAGS) $(TEST)/xAdiabaticPFR.cpp $(SOURCES)/AdiabaticPFR.o \
	$(LIBRARIES) -o $(BINS)/xAdiabaticPFR

xHeatWallPFR: $(OBJECTS)
	$(CXX) $(CPPFLAGS) $(TEST)/xHeatWallPFR.cpp $(SOURCES)/HeatWallPFR.o \
	$(LIBRARIES) -o $(BINS)/xHeatWallPFR

xIsothermalPFR: $(OBJECTS)
	$(CXX) $(CPPFLAGS) $(TEST)/xIsothermalPFR.cpp $(SOURCES)/IsothermalPFR.o \
	$(LIBRARIES) -o $(BINS)/xIsothermalPFR

clean:
	$(RM) $(OBJECTS) build

dist-clean: clean
	$(RM) $(BINS)/* CanteraPFR/lib/*
	cd CanteraPFR/doc && make clean

clean-all: dist-clean
	$(RM) external cygwin setup-x86_64.exe
