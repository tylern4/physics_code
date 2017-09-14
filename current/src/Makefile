UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
    FOPENMP = -fopenmp -lgfortran
endif

ROOTLIBS	= $(shell root-config --libs)
CXX = g++
CXXFLAGS = -O3 -march=native -fPIC -w -g $(FOPENMP) $(shell root-config --cflags)

SRCDIR   = ../src
OBJDIR   = ../obj
LIBDIR   = ../sobj

SRC = $(wildcard $(SRCDIR)/*.cpp)

OBJECTS  = $(SRC:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
LIBOUT =	$(LIBDIR)/lib.so

.PHONY: all clean

all: lib

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp
		$(CXX) $(CXXFLAGS) -c $< -o $@

sobj:
	@mkdir -p $(LIBDIR)

obj:
	@mkdir -p $(OBJDIR)

lib:	sobj $(OBJECTS)
	$(CXX) $(CXXFLAGS) -shared $(OBJECTS) $(ROOTLIBS) -o $(LIBOUT)

clean:
	-rm -f $(OBJECTS) $(LIBOUT)
