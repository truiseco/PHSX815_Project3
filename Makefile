# C++ compiler
CXX = c++

# necessary compiler flags for using ROOT (root.cern.ch) - remove these if you're not using root
ROOTCFLAGS    = $(shell root-config --cflags)
ROOTGLIBS     = $(shell root-config --glibs)

# ROOT shared library flags
#GLIBS = $(filter-out -stdlib=libc++ -pthread , $(ROOTGLIBS))
GLIBS = $(ROOTGLIBS)

# some compiler flags
CXXFLAGS = -std=c++17 -g -Wall -O2
# ROOT flags
CXXFLAGS += -fPIC $(filter-out -stdlib=libc++ -pthread , $(ROOTCFLAGS))
#CXXFLAGS = $(ROOTCFLAGS)

# location of source code
SRCDIR = ./src/

#location of header files
INCLUDEDIR = ./include/

CXXFLAGS += -I$(INCLUDEDIR)

# location of object files (from compiled library files)
OUTOBJ = ./obj/

CPP_FILES := $(wildcard src/*.cpp)
H_FILES := $(wildcard include/*.h)
OBJ_FILES := $(addprefix $(OUTOBJ),$(notdir $(CPP_FILES:.cpp=.o)))

# targets to make
all: ExpHypoSim.x ExpHypoTest.x

# recipe for building ExpHypoSim.x
ExpHypoSim.x:  $(SRCDIR)ExpHypoSim.C $(OBJ_FILES) $(H_FILES)
	$(CXX) $(CXXFLAGS) $(OUTOBJ)/*.o $ $< -o ExpHypoSim.x
	touch ExpHypoSim.x

# recipe for building ExpHypoTest.x
ExpHypoTest.x:  $(SRCDIR)ExpHypoTest.C $(OBJ_FILES) $(H_FILES)
	$(CXX) $(CXXFLAGS) $(OUTOBJ)/*.o $ $< $(GLIBS) -o ExpHypoTest.x
	touch ExpHypoTest.x

$(OUTOBJ)%.o: src/%.cpp include/%.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

# clean-up target (make clean)
clean:
	rm -f *.x
	rm -rf *.dSYM
	rm -f $(OUTOBJ)*.o
