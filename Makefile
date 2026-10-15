###############################################################################
##  CanteraPFR extension library Makefile
###############################################################################

CXX          := g++
AR           := ar
ARFLAGS      := rcs
RM           := rm -rf
MKDIR_P      := mkdir -p

OPTIONS      := -O3 -std=c++17 -Wall -Wextra
INCLUDES     := -Isrc/extension \
                -I/usr/local/include \
                -I/usr/include/eigen3 \
                -I/usr/lib/x86_64-linux-gnu/openmpi/include

LIBRARIES    := -L/usr/local/lib -Wl,-rpath,/usr/local/lib \
                -L/usr/lib/x86_64-linux-gnu/openmpi/lib \
                -lcantera_shared \
                -lsundials_ida \
                -lsundials_nvecserial \
                -lopenblas \
                -llapack \
                -lpthread

CPPFLAGS     := $(OPTIONS) -fPIC $(INCLUDES)

SRC_DIR      := src/extension
BUILD_DIR    := build
LIB_DIR      := lib

CPP_SRCS     := $(wildcard $(SRC_DIR)/*.cpp)
CPP_OBJS     := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(CPP_SRCS))

C_SRCS       := $(wildcard $(SRC_DIR)/*.c)
C_OBJS       := $(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(C_SRCS))

ALL_OBJS     := $(CPP_OBJS) $(C_OBJS)

STATIC_LIB   := $(LIB_DIR)/libCanteraPFR.a
SHARED_LIB   := $(LIB_DIR)/libCanteraPFR_shared.so

###############################################################################
##  TARGETS
###############################################################################

.PHONY: all static shared clean distclean

all: shared static

shared: $(SHARED_LIB)

static: $(STATIC_LIB)

$(STATIC_LIB): $(ALL_OBJS) | $(LIB_DIR)
	$(AR) $(ARFLAGS) $@ $^

$(SHARED_LIB): $(ALL_OBJS) | $(LIB_DIR)
	$(CXX) -shared -fPIC -o $@ $^ $(LIBRARIES)
	@$(MKDIR_P) CanteraPFR/lib
	@ln -sf ../../$@ CanteraPFR/lib/libCanteraPFR_shared.so

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) -x c++ -c $< -o $@

$(BUILD_DIR) $(LIB_DIR):
	$(MKDIR_P) $@

clean:
	$(RM) $(BUILD_DIR)

distclean: clean
	$(RM) $(LIB_DIR) CanteraPFR/lib

###############################################################################
##  EOF
###############################################################################
