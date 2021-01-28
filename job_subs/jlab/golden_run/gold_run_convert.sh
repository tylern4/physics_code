#!/bin/bash

#INPUT Sequences file

#SBATCH -J h10maker
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=256
#SBATCH -N 1
#SBATCH --output=log/h10maker_%A.out
#SBATCH --error=log/h10maker_%A.err


export NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUM_THREADS=2

export INPUT_FOLDER=/volatile/clas/clase1/tylern/bos_files
export OUTPUT_FOLDER=/work/clas/clase1/tylern/e1d/data/golden_run

mkdir -p $OUTPUT_FOLDER

echo $HOSTNAME

lscpu

echo "STARTING RUN"

START=`date +%s`

cat golden_run.list | parallel --dry-run -j$NUM_THREADS singularity exec /cvmfs/singularity.opensciencegrid.org/tylern4/clas6:latest h10maker -r $INPUT_FOLDER/{} $OUTPUT_FOLDER/h10_gr_{.}.root

END=`date +%s`
echo Execution time was `expr $end - $start` seconds.


echo "DONE!!"
