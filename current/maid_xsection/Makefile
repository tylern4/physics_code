PROGRAM=example
SOURCE=example.cc

# Compilation and linking settings:
FCFLAGS = -O3 -m64 -fno-automatic -ffixed-line-length-none -fno-second-underscore
LIB_PATH = -L$(CERNLIB) -L$(CLAS6)/lib
LIB = -lmaid_xsection -lgfortran -lrecutl -lmathlib -lkernlib $(shell root-config --libs)
CXXFLAGS = $(shell root-config --cflags)
CXX = $(shell root-config --cxx)

# Targets:
$(PROGRAM): $(PROGRAM).o
	$(CXX) $(CXXFLAGS) $(PROGRAM).o $(LIB_PATH) $(LIB) -o $(PROGRAM)

$(PROGRAM).o: $(SOURCE)
	$(CXX) $(CXXFLAGS) -c $(SOURCE) -o $(PROGRAM).o

clean:
	$(RM) $(PROGRAM) $(PROGRAM).o

