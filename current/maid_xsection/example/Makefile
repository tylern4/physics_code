# This Makefile shows how to compile a program called "example" which
# consists of a single source code file, example.cc.  Simply change
# the variables PROGRAM and SOURCE to match what you want for your
# program.
#
# The key pieces to note are
#
# 1. gfortran must be used to link the code.  I've tried to find ways
#    to make g++ link to the Fortran libraries appropriately, but it
#    is more difficult than it's worth.
#
# 2. Keep the order of the gfortran flags, library location flags,
#    etc. as they are in the program compilation recipe.  I have no
#    idea why, but changing the order sometimes breaks the compilation
#    when I tested it.  This is the order I've found to work in
#    general.

# Program name and source code path
PROGRAM=example
SOURCE=example.cc

# Compilation and linking settings:
FORTRAN_FLAGS = -O3 -m64 -fno-automatic -ffixed-line-length-none -fno-second-underscore
#LIBRARYLOC_FLAGS= -L$(HOME)/slib/Linux64RHEL5 -L$(HOME)/lib/Linux64RHEL5 -L/group/clas/builds/PRODUCTION/slib/Linux64RHEL5 -L/u/home/gleb/lib/LinuxRHFC8 -L/apps/cernlib/x86_64_rhel7/2005/lib -L/usr/lib -L/usr/X11R6/lib64/ -L/usr/lib/mysql
#LIBRARY_FLAGS = -lclasutil -lfputil -lrecutl -lbos -lfpack -lonline_dummy -lc_bos_io -llepto651 -ljetset74 -lpacklib -lmathlib -lkernlib -llepto651 -lfputil -lbos -lfpack -lnsl -ldl -lg2c $(shell root-config --libs)
COMPILE_FLAGS = $(shell root-config --cflags)

# Targets:
$(PROGRAM): $(PROGRAM).o
	g++  $(COMPILE_FLAGS) $(shell /usr/bin/cernlib) $(PROGRAM).o  -L.. -lmaid_xsection $(shell root-config --libs) -lminicern -o $(PROGRAM)

$(PROGRAM).o: $(SOURCE)
	g++ $(shell /usr/bin/cernlib) $(COMPILE_FLAGS) -c $(SOURCE) -o $(PROGRAM).o

clean:
	rm -f $(PROGRAM) $(PROGRAM).o
