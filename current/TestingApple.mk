ROOTLIBS	= $(shell root-config --libs)
INCLUDE		= $(shell root-config --incdir)
ROOTSYS		= $(shell root-config --exec-prefix)
CXX = clang++-git
#CXX = /Users/tylern/Downloads/usr/local/bin/g++
CXXFLAGS =      -O2 -fPIC -w -g $(shell root-config --cflags) -frtti
INCS =          -I$(shell root-config --incdir)
LIBS =          $(shell root-config --glibs)
TARGET =	MacE1d

all:	clean E1D

E1D:	main.o 
	$(CXX) main.o -L. $(CXXFLAGS) $(LIBS) -o $(TARGET) 

clean:
	-rm -f $(OBJS) $(LIBOUT) MacE1d main.o 
