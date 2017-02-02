UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
    FOPENMP = -fopenmp -lgfortran
endif

ROOTLIBS	= $(shell root-config --libs)
CXX = g++
CXXFLAGS =      -O2 -fPIC -w -g $(FOPENMP) $(shell root-config --cflags)
TARGET =	    e1d
SRCDIR   = ../src
OBJDIR   = ../obj
LIBDIR   = ../sobj
BINDIR   = ../bin

SRC = $(wildcard $(SRCDIR)/*.cpp) $(wildcard *.cpp)
OBJ=	$(patsubst %.cpp,%.o,$(SRC))
LIBOUT =	$(LIBDIR)/lib.so

.PHONY: all clean

all:	clean $(BINDIR)/$(TARGET)

#$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
#	$(CXX) $(CXXFLAGS)  -c $< -o $@

$(BINDIR)/$(TARGET): obj $(OBJ)
	$(CXX) $(OBJ) $(CXXFLAGS) $(ROOTLIBS) -o $(TARGET)

sobj:
	@mkdir -p $@

obj:
	@mkdir -p $@

lib:	sobj $(OBJ)
	$(CXX) $(CXXFLAGS) -shared $(OBJ) $(ROOTLIBS) -o $(LIBOUT)

clean:
	-rm -f $(TARGET) $(OBJ)
