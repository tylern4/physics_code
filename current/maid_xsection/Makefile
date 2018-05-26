FC = gfortran
CXX = $(shell root-config --cxx)
UNAME := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
	LDFLAGS = -lstdc++
else
	LDFLAGS = -lc++
endif

ROOTLIBS	= $(shell root-config --libs)
CXXFLAGS = -O3 -Wall -Wextra -march=native -fPIC -w -g $(shell root-config --cflags)
FCFLAGS = -std=legacy -fPIC -I.

C_SRC = $(wildcard *.cc)
F_SRC = maid.F splie2.F maid_xsection.F spline.F splin2.F splint.F dsigma.F aaosig1.F aao.F daresbury.F maid_lee.F 
#$(wildcard *.F) helicity_amps.F read_sf_file.F cgln_amps.F legendre.F multipole_amps.F xsection.F interp.F
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
	#$(CXX) $(CXXFLAGS) -shared $(C_OBJECTS) $(F_OBJECTS) $(ROOTLIBS) $(LDFLAGS) -o $(LIBOUT)


clean :
	$(RM) -rf *.o libmaid_xsection.so
