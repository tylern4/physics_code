#!/bin/bash

#INPUT Sequences file

#SBATCH -J clas12_analysis
#SBATCH -n 20
#SBATCH -N 1
#SBATCH -p msmoms
#SBATCH --output=log/clas6_%a.out
#SBATCH --error=log/clas6_%a.err

module load gcc/6.4.0
module load python3/anaconda/5.2.0

export CC=$(which gcc)
export CXX=$(which g++)

#CMAKE
export CMAKE_DIR=/work/gothelab/software/cmake
export PATH=$CMAKE_DIR/bin:$PATH

#CERN ROOT
export ROOTSYS=/work/gothelab/software/root
export PATH=$ROOTSYS/bin:$PATH
export PYTHONDIR=$ROOTSYS
export LD_LIBRARY_PATH=$ROOTSYS/lib:$PYTHONDIR/lib:$ROOTSYS/bindings/pyroot:$LD_LIBRARY_PATH
export PYTHONPATH=/usr/local/lib:$ROOTSYS/lib:$PYTHONPATH:$ROOTSYS/bindings/pyroot

cd /home/tylerns/physics_code/current/build

./e1d /work/gothelab/clas6/e1d/e1d_v2.root /work/gothelab/clas6/e1d/data/v2/h10_r2*

./e1d_mc /work/gothelab/clas6/e1d/e1d_sim.root /work/gothelab/clas6/e1d/sim/root/sim*.root
