UNAME := $(shell uname)
ifeq ($(UNAME), Linux)
    FOPENMP = -fopenmp -lgfortran
endif

ROOTLIBS	= $(shell root-config --libs)
CXX = g++ -stdlib=libc++
CXXFLAGS =      -O2 -fPIC -w -g $(FOPENMP) $(shell root-config --cflags) 
TARGET =	    e1d
SRC =		$(wildcard *.cpp)
OBJ=	$(patsubst %.cpp,obj/%.o,$(SRC))  
#-lgfortran
LIBOUT =        sobj/lib.so


.PHONY: all clean

all:	$(TARGET)

obj/%.o: %.cpp
	$(CXX) $(CXXFLAGS)  -c $< -o $@

$(TARGET): obj $(OBJ) 
	$(CXX) $(OBJ) -L. $(CXXFLAGS) $(ROOTLIBS) -o $(TARGET) 

sobj:
	@mkdir -p $@

obj:
	@mkdir -p $@

lib:	sobj $(OBJ)
	$(CXX) $(CXXFLAGS) -shared $(OBJ) $(ROOTLIBS) -o $(LIBOUT)

clean:
	-rm -rf $(TARGET) $(OBJ) obj/
