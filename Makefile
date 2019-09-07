# Define any paths that are external to the current directory here, this boost path works only for lxplus 
# account holders (CERN members) but can be easily modified to your local installation of boost
BOOST_PATH = /usr/local/include/boost_1_62_0/
# Input the title of your "main" file here, e.g. test.cpp.
main := main

CXXFLAGS	:= -Wall -fPIC -W -g -std=c++11 -pedantic

BOOSTLIBS	:= -L $(BOOST_PATH)/lib
CXXLIBS		:= -L.

CXX	:= g++ -g -L/usr/local/boost_1_62_0/stage/lib/ -lboost_filesystem -lboost_system

# We could place here different extensions such as *.hh, *.h, *.H etc
INCLUDES := -I./include/$(wildcard *.hpp)		# This would contain any -I options to the
INCLUDES += -I$(BOOST_PATH)/include			# compiler, if there are any e.g. Boost and local

# Tell it the location of the src directory
SRC_MODULES   := src
SRC_DIR   := $(addprefix ${PWD}/,$(SRC_MODULES))
# Tell it to look also in current directory for the <"main">.cpp file
SRC_DIR	  += ${PWD}

BUILD_DIR := bin

# Add all locations of source files. We could place here different extensions such as *.cc, *.c, *.C etc
SRC       := $(foreach sdir,$(SRC_DIR),$(wildcard $(sdir)/*.cpp))
# Create associated objects for each in some 
OBJ       := $(patsubst $(SRC)/%.cpp,bin/%.o,$(SRC))

# Again could change extensions
vpath %.cpp $(SRC_DIR)


.PHONY: all checkdirs clean

all: checkdirs bin/$(main).exe


bin/$(main).exe: $(OBJ)
	@echo "Linking Exectuables ..."
	$(CXX) $(CXXFLAGS) $(INCLUDES) $^ -o $@ $(CXXLIBS)
	@echo "Done."


checkdirs: $(BUILD_DIR)

$(BUILD_DIR):
	@mkdir -p $@

# Clean up ..
clean:
	@echo "Cleaning up...."
	@rm -rf $(BUILD_DIR) 
	@echo "Done."



