#!/bin/bash

#INPUT Sequences file

#SBATCH -J e1d_csv_maker
#SBATCH -n 40
#SBATCH -N 1
#SBATCH -p BigMem
#SBATCH --output=log/csv_maker_%a_%A.out
#SBATCH --error=log/csv_maker_%a_%A.err

module load gcc/6.4.0
module load python3/anaconda/5.2.0

export CC=$(which gcc)
export CXX=$(which g++)


#CERN ROOT
export ROOTSYS=/work/gothelab/software/root
export PATH=$ROOTSYS/bin:$PATH
export PYTHONDIR=$ROOTSYS
export LD_LIBRARY_PATH=$ROOTSYS/lib:$PYTHONDIR/lib:$ROOTSYS/bindings/pyroot:$LD_LIBRARY_PATH
export PYTHONPATH=/usr/local/lib:$ROOTSYS/lib:$PYTHONPATH:$ROOTSYS/bindings/pyroot

cd /home/tylerns/physics_code/build

export NUM_THREADS=40

./csv_maker_mc -o /work/tylerns/e1d/outputs/csv/mc_old -e1d /work/gothelab/clas6/e1d/sim/osg
cat /work/tylerns/e1d/outputs/csv/mc_old_*_e1d.csv > /work/tylerns/e1d/outputs/csv/mc_old_e1d.csv
xz -3 -v -T4 /work/tylerns/e1d/outputs/csv/mc_old_e1d.csv
rm -rf /work/tylerns/e1d/outputs/csv/mc_old_*_e1d.csv

