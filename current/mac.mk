ROOTLIBS	= $(shell root-config --libs)
INCLUDE		= $(shell root-config --incdir)
ROOTSYS		= $(shell root-config --exec-prefix)
C++ = clang++
CXX = clang++
CXXFLAGS =      -O2 -fPIC -w -g -fopenmp $(shell root-config --cflags) 
#CXXFLAGS =		-O2 -fPIC -w -g -fopenmp -pthread -m64 -I/usr/local/root/include
INCS =          -I$(shell root-config --incdir)
LIBS =          $(shell root-config --glibs)
TARGET =	    e1d

all:	clean E1D

E1D:	main.o 
	$(CXX) main.o -L. $(CXXFLAGS) $(LIBS) -o $(TARGET) 

clean:
	-rm -f $(OBJS) $(LIBOUT) e1d main.o 

#lib:	$(OBJS)
#	$(CXX) $(CXXFLAGS) -shared $(OBJS) 
