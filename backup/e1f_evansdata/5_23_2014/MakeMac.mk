ROOTLIBS	= $(shell root-config --libs)
INCLUDE		= $(shell root-config --incdir)
ROOTSYS		= $(shell root-config --exec-prefix)
C++ = clang
CXX = clang++
CXXFLAGS =      -O2 -fPIC -w -g -fopenmp $(shell root-config --cflags) 
INCS =          -I$(shell root-config --incdir)
LIBS =          $(shell root-config --glibs)
TARGET =	    test

all:	clean lib $(TARGET)

%.o: %.cpp
	$(CXX) $(INCS) $(CXXFLAGS) -c test.cpp

$(TARGET):	lib test.o
	$(CXX) test.o -L. $(CXXFLAGS) $(LIBS) -o $(TARGET) 

clean:
	-rm -f $(OBJS) $(LIBOUT) test test.o a.out

lib:	$(OBJS)
	$(CXX) $(CXXFLAGS) -shared $(OBJS) $(LIBS) 
