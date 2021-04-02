#!/bin/bash

#INPUT Sequences file

#SBATCH -J e1d_mc
#SBATCH -n 28
#SBATCH -N 1
#SBATCH -p defq,defq-48core,BigMem,msmoms
#SBATCH --output=log/e1d_mc_%A.out
#SBATCH --error=log/e1d_mc_%A.err

module load gcc/6.4.0
module load python3/anaconda/2020.02

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

cd /home/tylerns/physics_code/build

export NUM_THREADS=$SLURM_JOB_CPUS_PER_NODE 
echo $SLURM_JOB_CPUS_PER_NODE

./e1d_mc -o /work/tylerns/e1d/outputs/e1d_sim.root -i /work/gothelab/clas6/e1d/sim/npip_test

