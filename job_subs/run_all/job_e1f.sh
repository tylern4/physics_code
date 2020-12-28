#!/bin/bash

#INPUT Sequences file

#SBATCH -J e1f_data
#SBATCH -n 28
#SBATCH -N 1
#SBATCH -p defq,BigMem,defq-48core,gpu-v100-16gb,gpu-v100-32gb,msmoms
#SBATCH --output=log/e1f_data_%A.out
#SBATCH --error=log/e1f_data_%A.err

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

./e1f /work/tylerns/e1f/outputs/e1f_new.root /work/gothelab/clas6/e1f/skimmed/pass2/v1/*.root

./csv_maker -o /work/tylerns/e1f/outputs/csv/data -e1f /work/gothelab/clas6/e1f/skimmed/pass2/v1

cat /work/tylerns/e1f/outputs/csv/data_*_e1f.csv > /work/tylerns/e1f/outputs/csv/data_e1f.csv
rm -rf /work/tylerns/e1f/outputs/csv/data_*_e1f.csv
/work/apps/python3/anaconda3/2020.02/bin/python /home/tylerns/physics_code/python/csv_to_feather.py /work/tylerns/e1f/outputs/csv/data_e1f.csv

