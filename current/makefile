ROOTLIBS	= $(shell root-config --libs)
INCLUDE		= $(shell root-config --incdir)
ROOTSYS		= $(shell root-config --exec-prefix)
CXX = g++
CXXFLAGS =      -O2 -fPIC -w -g $(shell root-config --cflags) 
INCS =          -I$(shell root-config --incdir)
LIBS =          $(shell root-config --glibs)
TARGET =	    e1d
SRC =		$(wildcard *.cpp)
FILENAME=	$(SRC:.cpp=.o)

all:	clean E1D

E1D: $(FILENAME) 
	$(CXX) $(FILENAME) -L. $(CXXFLAGS) $(LIBS) -o $(TARGET) 

clean:
	-rm -f $(OBJS) $(LIBOUT) $(TARGET) $(FILENAME)
