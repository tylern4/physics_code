#!/bin/bash

#INPUT Sequences file

#SBATCH -J e1d_csv_maker
#SBATCH -n 16
#SBATCH -N 1
#SBATCH -p defq-48core,gpu-v100-16gb
##SBATCH -p defq-48core,gpu-v100-16gb,gpu-v100-32gb,msmoms,defq,gpu
##SBATCH -p msmoms
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

export NUM_THREADS=$SLURM_JOB_CPUS_PER_NODE


if [ "$HOSTNAME" = ultra ]; then
    echo "Using ultra settings"
    export OUTPUT_DIR=/data3/scratch/$USER/mc
else
    export OUTPUT_DIR=/local/$USER/mc
fi

mkdir -p $OUTPUT_DIR

echo $NUM_THREADS
echo $HOSTNAME


res1=$(date +%s.%N)

./csv_maker_mc -o $OUTPUT_DIR/mc_osg -e1d /work/gothelab/clas6/e1d/sim/npip

mv $OUTPUT_DIR/mc_osg*.csv /work/tylerns/e1d/outputs/csv
rm -rf $OUTPUT_DIR/*

res2=$(date +%s.%N)
dt=$(echo "$res2 - $res1" | bc)
dd=$(echo "$dt/86400" | bc)
dt2=$(echo "$dt-86400*$dd" | bc)
dh=$(echo "$dt2/3600" | bc)
dt3=$(echo "$dt2-3600*$dh" | bc)
dm=$(echo "$dt3/60" | bc)
ds=$(echo "$dt3-60*$dm" | bc)
printf "Total runtime copy: %d:%02d:%02d:%02.4f\n" $dd $dh $dm $ds


#cat /work/tylerns/e1d/outputs/csv/mc_osg_*_e1d.csv > /work/tylerns/e1d/outputs/csv/mc_osg_e1d.csv
#cp /work/tylerns/e1d/outputs/csv/mc_osg_0_e1d.csv /work/tylerns/e1d/outputs/csv/mc_test_e1d.csv
#rm -rf /work/tylerns/e1d/outputs/csv/mc_osg_*_e1d.csv

