# Makefile for preplot
# Dependencies generated automatically.
#
# Default compiler is icpc
CXX := icpc
#CXX := /usr/local/Cellar/gcc/5.3.0/bin/g++-5
DEBUG = -g
CXXFLAGS = -gcc-name=/usr/bin/gcc -openmp -std=c++11
#CXXFLAGS = -fopenmp -std=c++11

# Compile all cpp files in source directory
SRCS := $(shell find . -name "*.cpp")

# Substitute source files for object files
OBJS := $(subst .o,.cpp,$(SRCS))

# default target
all : preplot

preplot : $(OBJS)
	$(CXX) $(CXXFLAGS) $(DEBUG) -o preplot $(OBJS)

depend: .depend

.depend: $(SRCS)
	rm -f ./.depend
	$(CXX) $(CXXFLAGS) -MM $^>>./.depend;

clean:
	rm -f $(OBJS)

include .depend
