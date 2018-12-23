###############################################################################
##  MAKEFILE
###############################################################################

SOURCES   := src
SUNDIALS  := -lsundials_ida -lsundials_nvecserial
INCLUDES  := -Iexternal/include
LIBRARIES := -Lexternal/lib -lcantera_shared -pthread $(SUNDIALS)
CPPFLAGS  := -O3 -std=c++11 $(INCLUDES)
CXX       := g++

###############################################################################
##  CONVERT OBJECTS
###############################################################################

SOBJECT := $(shell find $(SOURCES) -iname "*.cpp")
OBJECTS := $(patsubst $(SOURCES)/%.cpp,$(SOURCES)/%.o,$(SOBJECT))

###############################################################################
##  COMPILE MODELS
## TODO for now just the available model, later this has to be improved!
###############################################################################

all: $(OBJECTS)
	$(CXX) $(CPPFLAGS) $(OBJECTS) $(LIBRARIES) -o IsothermalPFR
