UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
    FOPENMP = -fopenmp
endif

ROOTLIBS	= $(shell root-config --libs)
CXX = g++
CXXFLAGS = -O3 -march=native -fPIC -w -g $(FOPENMP) $(shell root-config --cflags)
TARGET =	    511_lab
SRCDIR   = ../src
OBJDIR   = ../obj
LIBDIR   = ../sobj
BINDIR   = ../bin

SRC = $(wildcard $(SRCDIR)/*.cpp)
MAINS = ./511_lab.cpp

OBJECTS  = $(SRC:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
MAIN = $(MAINS:./%.cpp=$(OBJDIR)/%.o)
LIBOUT =	$(LIBDIR)/lib.so

.PHONY: all clean

all:	$(BINDIR)/$(TARGET)

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp
		$(CXX) $(CXXFLAGS) -c $< -o $@

$(MAIN): $(OBJDIR)/%.o : ./%.cpp
		$(CXX) $(CXXFLAGS) -c $< -o $@

$(BINDIR)/$(TARGET): obj $(OBJECTS) $(MAIN)
	$(CXX) $(OBJECTS) $(MAIN) $(CXXFLAGS) $(ROOTLIBS) -o $(BINDIR)/$(TARGET)

sobj:
	@mkdir -p $(LIBDIR)

obj:
	@mkdir -p $(OBJDIR)

lib:	sobj $(OBJECTS)
	$(CXX) $(CXXFLAGS) -shared $(OBJECTS) $(ROOTLIBS) -o $(LIBOUT)

semiclean:
	-rm -f $(MAIN)

clean:
	-rm -f $(BINDIR)/$(TARGET) $(OBJECTS) $(MAIN)
