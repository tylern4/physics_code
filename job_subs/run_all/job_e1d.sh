#!/bin/bash

#INPUT Sequences file

#SBATCH -J e1d_data
#SBATCH -n 28
#SBATCH -N 1
#SBATCH -p defq,BigMem,gpu,defq-48core,gpu-v100-16gb,gpu-v100-32gb,msmoms
#SBATCH --output=log/e1d_data_%A.out
#SBATCH --error=log/e1d_data_%A.err

module load gcc/6.4.0
module load python3/anaconda/2020.02
module load valgrind/3.14

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

#valgrind --tool=callgrind --callgrind-out-file=/work/tylerns/e1d_$SLURM_JOB_ID.callgrind.out ./e1d /work/tylerns/e1d/outputs/e1d_new.root /work/gothelab/clas6/e1d/skimmed/pass2/v2/full/*.root

# ./e1d /work/tylerns/e1d/outputs/e1d_skim.root /work/gothelab/clas6/e1d/skimmed/pass2/v2/full/*.root
./e1d /work/tylerns/e1d/outputs/e1d_data.root /work/gothelab/clas6/e1d/golden_run/*.root
./e1d /work/tylerns/e1d/outputs/e1d_empty.root /work/gothelab/clas6/e1d/empty/*.root

# ./csv_maker -o /work/tylerns/e1d/outputs/csv/skim -e1d /work/gothelab/clas6/e1d/skimmed/pass2/v2/full
./csv_maker -o /work/tylerns/e1d/outputs/csv/data -e1d /work/gothelab/clas6/e1d/golden_run
./csv_maker -o /work/tylerns/e1d/outputs/csv/empty -e1d /work/gothelab/clas6/e1d/empty


#cat /work/tylerns/e1d/outputs/csv/data_*_e1d.csv > /work/tylerns/e1d/outputs/csv/data_e1d.csv
#rm -rf /work/tylerns/e1d/outputs/csv/data_*_e1d.csv
#/work/apps/python3/anaconda3/2020.02/bin/python /home/tylerns/physics_code/python/csv_to_feather.py /work/tylerns/e1d/outputs/csv/data_e1d.csv

