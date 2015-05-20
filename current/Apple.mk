ROOTLIBS	= $(shell root-config --libs)
INCLUDE		= $(shell root-config --incdir)
ROOTSYS		= $(shell root-config --exec-prefix)
CXX = clang++
CXXFLAGS =      -O2 -fPIC -w -g $(shell root-config --cflags) 
INCS =          -I$(shell root-config --incdir)
LIBS =          $(shell root-config --glibs)
TARGET =		MacE1d
FILENAME=		main

all:	clean E1D

E1D:	$(FILENAME).o 
	$(CXX) $(FILENAME).o -L. $(CXXFLAGS) $(LIBS) -o $(TARGET) 

clean:
	-rm -f $(OBJS) $(LIBOUT) $(TARGET) $(FILENAME).o 
