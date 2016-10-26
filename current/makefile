UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
FOPENMP = -fopenmp
endif

ROOTLIBS	= $(shell root-config --libs)

CXX = g++
CXXFLAGS =      -O2 -fPIC -w -g $(FOPENMP) $(shell root-config --cflags) 
TARGET =	    e1d
SRC =		$(wildcard *.cpp)
OBJ=	$(SRC:.cpp=.o)

.PHONY: all clean

all:	clean $(TARGET)

$(TARGET): $(OBJ) 
	$(CXX) $(OBJ) -L. $(CXXFLAGS) $(ROOTLIBS) -o $(TARGET) 

clean:
	-rm -f $(TARGET) $(OBJ)
