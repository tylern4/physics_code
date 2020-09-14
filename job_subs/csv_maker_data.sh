#!/bin/bash

#INPUT Sequences file

#SBATCH -J e1d_csv_maker
#SBATCH -n 48
#SBATCH -N 1
#SBATCH -p defq-48core,gpu-v100-16gb,gpu-v100-32gb
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

export NUM_THREADS=$SLURM_JOB_CPUS_PER_NODE


echo $HOSTNAME

./csv_maker -o /work/tylerns/e1d/outputs/csv/data -e1d /work/gothelab/clas6/e1d/data/v2 -e1f /work/gothelab/clas6/e1f/data
cat /work/tylerns/e1d/outputs/csv/data_*_e1d.csv > /work/tylerns/e1d/outputs/csv/data_e1d.csv
rm -rf /work/tylerns/e1d/outputs/csv/data_*_e1d.csv

/work/apps/python3/anaconda3/2020.02/bin/python /home/tylerns/physics_code/python/csv_to_feather.py /work/tylerns/e1d/outputs/csv/data_e1d.csv

#cat /work/tylerns/e1d/outputs/csv/data_*_e1f.csv > /work/tylerns/e1d/outputs/csv/data_e1f.csv
#rm -rf /work/tylerns/e1d/outputs/csv/data_*_e1f.csv

#/work/apps/python3/anaconda3/2020.02/bin/python /home/tylerns/physics_code/python/csv_to_feather.py /work/tylerns/e1d/outputs/csv/data_e1f.csv
