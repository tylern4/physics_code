UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
    FOPENMP = -fopenmp 
endif
ROOTLIBS	= $(shell root-config --libs) -lCling
CXX = $(shell root-config --cxx)
CXXFLAGS = -O3 -rdynamic -shared -march=native -fPIC -w -g $(FOPENMP) $(shell root-config --cflags)
INC = $(root-config --incdir)

SRCDIR   = ../src
OBJDIR   = ../obj

SRCS = $(wildcard $(SRCDIR)/*.cpp)
OBJECTS  = $(SRCS:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)

GEN = genreflex
SRC = $(wildcard *.h)
PYC = $(wildcard *.pyc)
CPP = $(SRC:%.h=%_rflx.cpp)
SOBJ  = $(SRC:%.h=%.so)
TEMPS = $(wildcard *.pcm) $(wildcard *_tmp) $(wildcard *.so.dSYM) $(wildcard *_rflx*)

.PHONY: all clean

all: $(OBJECTS) $(CPP) $(SOBJ)

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(CPP): $(SRC)
	$(foreach s,$(SRC),$(GEN) $(s);)

$(SOBJ): %.so : %_rflx.cpp
		$(CXX) $(OBJECTS) $(CXXFLAGS) -I$(INC) -c $< -o $@ $(ROOTLIBS)


clean:
	-rm -rf $(SOBJ) $(OBJECTS) $(TEMPS) $(PYC)
