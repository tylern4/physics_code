ROOTLIBS	= $(shell root-config --libs)
INCLUDE		= $(shell root-config --incdir)
ROOTSYS		= $(shell root-config --exec-prefix)
C++ = g++
CXX = g++
CXXFLAGS =      -O2 -fPIC -w -g -fopenmp $(shell root-config --cflags) 
CLAS_INC=		/home/tylern/clas/ClasTool/include
INCS =          -I$(shell root-config --incdir) -I$(CLAS_INC)
LIBS =          $(shell root-config --glibs)
TARGET =	    convertRoot

all:	clean convert

convert:	convertRoot.o 
	$(CXX) convertRoot.o -L. $(CXXFLAGS) $(LIBS) -o $(TARGET) 

clean:
	-rm -f $(OBJS) $(LIBOUT) convertRoot convertRoot.o 

lib:	$(OBJS)
	$(CXX) $(CXXFLAGS) -shared $(OBJS) $(LIBS) 
