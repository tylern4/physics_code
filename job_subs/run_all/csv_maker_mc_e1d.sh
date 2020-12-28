#!/bin/bash

#INPUT Sequences file

#SBATCH -J e1d_csv_maker
#SBATCH -n 48
#SBATCH -N 1
#SBATCH -p defq-48core,msmoms,gpu-v100-16gb,gpu-v100-32gb
#SBATCH --output=log/csv_maker_%a_%A.out
#SBATCH --error=log/csv_maker_%a_%A.err

module load gcc/6.4.0

export CC=$(which gcc)
export CXX=$(which g++)


#CERN ROOT
export ROOTSYS=/work/gothelab/software/root
export PATH=$ROOTSYS/bin:$PATH
export PYTHONDIR=$ROOTSYS
export LD_LIBRARY_PATH=$ROOTSYS/lib:$PYTHONDIR/lib:$ROOTSYS/bindings/pyroot:$LD_LIBRARY_PATH
export PYTHONPATH=/usr/local/lib:$ROOTSYS/lib:$PYTHONPATH:$ROOTSYS/bindings/pyroot


cd /home/tylerns/physics_code/build


#export NUM_THREADS=$(( $SLURM_JOB_CPUS_PER_NODE * 2))
export NUM_THREADS=$SLURM_JOB_CPUS_PER_NODE

echo $NUM_THREADS
echo $HOSTNAME

./csv_maker_mc -o /work/tylerns/e1d/outputs/csv/mc_osg -e1d /work/gothelab/clas6/e1d/sim/npip

#cat /work/tylerns/e1d/outputs/csv/mc_osg_*_e1d.csv > /work/tylerns/e1d/outputs/csv/mc_osg_e1d.csv
#cp /work/tylerns/e1d/outputs/csv/mc_osg_0_e1d.csv /work/tylerns/e1d/outputs/csv/mc_test_e1d.csv
#rm -rf /work/tylerns/e1d/outputs/csv/mc_osg_*_e1d.csv

