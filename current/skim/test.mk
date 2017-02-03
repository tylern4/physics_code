UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
    FOPENMP = -fopenmp -lgfortran
endif

ROOTLIBS	= $(shell root-config --libs)
CXX = g++
CXXFLAGS =      -O2 -fPIC -w -g $(FOPENMP) $(shell root-config --cflags)
TARGET = Skim
SRCDIR   = ../src
OBJDIR   = ../obj
LIBDIR   = ../sobj
BINDIR   = ../bin

SRC = $(wildcard $(SRCDIR)/*.cpp)
MAINS = $(wildcard ./*.cpp)

OBJECTS  = $(SRC:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
MAIN = $(MAINS:./%.cpp=$(OBJDIR)/%.o)
LIBOUT =	$(LIBDIR)/lib.so

.PHONY: all clean

all:	clean $(BINDIR)/$(TARGET)

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

clean:
	-rm -f $(BINDIR)/$(TARGET) $(OBJECTS) $(MAIN)
