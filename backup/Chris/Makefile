ROOTLIBS	= $(shell root-config --libs)
INCLUDE		= $(shell root-config --incdir)
ROOTSYS		= $(shell root-config --exec-prefix)
CXX = g++
CXXFLAGS =      -O2 -fPIC -w -g $(shell root-config --cflags) 
INCS =          -I$(shell root-config --incdir)
LIBS =          $(shell root-config --glibs)
#This is what your programs is going to be named
TARGET =	attempt8
#This is what your file is named without the .C 
FILENAME =	attempt8

all:	clean test

test:	$(FILENAME).o 
	$(CXX) $(FILENAME).o -L. $(CXXFLAGS) $(LIBS) -o $(TARGET) 

clean:
	-rm -f $(OBJS) $(LIBOUT) $(TARGET) $(FILENAME).o 
