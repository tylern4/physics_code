ROOTLIBS	= $(shell root-config --libs)
INCLUDE		= $(shell root-config --incdir)
ROOTSYS		= $(shell root-config --exec-prefix)
CXX = clang++
CXXFLAGS =      -O2 -fPIC -w -g $(shell root-config --cflags) -emit-llvm -S
INCS =          -I$(shell root-config --incdir)
LIBS =          $(shell root-config --glibs)
TARGET =	    Skim
SRC =		$(wildcard *.cpp)
FILENAME=	$(SRC:.cpp=.o)

all:	clean SKIM

SKIM:	$(FILENAME) 
	$(CXX) $(FILENAME) -L. $(CXXFLAGS) $(LIBS) -o $(TARGET) 

clean:
	-rm -f $(OBJS) $(LIBOUT) $(FILENAME)