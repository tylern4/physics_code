FC = gfortran
CXX = $(shell root-config --cxx)
UNAME := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
	LDFLAGS = -lstdc++
else
	LDFLAGS = -lc++
endif

CXXFLAGS = -O3 -fPIC
FCFLAGS = -std=legacy -fPIC -I.

C_SRC = $(wildcard *.cc)
F_SRC = $(wildcard *.F)

C_OBJECTS  = $(C_SRC:%.cc=%.o)
F_OBJECTS = $(F_SRC:%.F=%.o)

LIBOUT =	libmaid_xsection.so

.PHONY: all clean

all: lib

$(C_OBJECTS) : %.o : %.cc
		$(CXX) $(CXXFLAGS) -c $< -o $@

$(F_OBJECTS) : %.o : %.F
		$(FC) $(FCFLAGS) -c $< -o $@

lib:	$(C_OBJECTS) $(F_OBJECTS)
	$(FC) $(FCFLAGS) -shared $(C_OBJECTS) $(F_OBJECTS) $(LDFLAGS) -o $(LIBOUT)

clean :
	$(RM) -rf *.o libmaid_xsection.so
