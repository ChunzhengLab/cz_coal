# Makefile for ROOT-based program using g++

# Compiler settings
CXX := g++
CXXFLAGS := $(shell root-config --cflags) -O2 -Wall
LDFLAGS := $(shell root-config --libs --glibs) -Wl,-rpath,$(shell root-config --libdir)

# Executable name
EXE := randPartonGen

# Source files
SRC := RandPartonGenerator.cpp

# Object files
OBJ := $(SRC:.cpp=.o)

# Build rule for the executable
$(EXE): $(OBJ)
	$(CXX) -o $@ $^ $(LDFLAGS)

# Rule to compile object files
%.o: %.cpp
	$(CXX) -c $< $(CXXFLAGS)

# Clean up
.PHONY: clean
clean:
	rm -f $(OBJ) $(EXE)