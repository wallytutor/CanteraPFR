###############################################################################
##  MAKEFILE
###############################################################################

include Makefile.in

SUNDIALS  := -L$(SUNDIALS_PATH) $(SUNDIALS_LIBS)
CANTERA   := -L$(CANTERA_PATH) $(CANTERA_LIBS)
INCLUDES  := -I$(EXTERNAL)/include -ICanteraPFR/include
LIBRARIES := $(CANTERA) $(SUNDIALS) $(OPENBLAS) -pthread
CPPFLAGS  := $(OPTIONS) -fPIC $(INCLUDES)

CXX       := g++
AR        := ar
ARFLAGS   := sq
RM        := rm -rf

SOURCES   := CanteraPFR/src/cpp
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
	python3 setup.py install && cd CanteraPFR/doc; make html

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
