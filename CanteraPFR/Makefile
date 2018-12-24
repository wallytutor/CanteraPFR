###############################################################################
##  MAKEFILE
###############################################################################

include Makefile.in

SOURCES   := src
SUNDIALS  := -lsundials_ida -lsundials_nvecserial
INCLUDES  := -Iinclude -I$(EXTERNAL)/include
# FIXME shared library failed on Cygwin!
LIBRARIES := -L$(EXTERNAL)/lib -lcantera -pthread $(SUNDIALS)
CPPFLAGS  := $(OPTIONS) $(INCLUDES)
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
