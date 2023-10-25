#!/bin/bash

#INPUT Sequences file

#SBATCH -J h10maker
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=256
#SBATCH -N 1
#SBATCH --output=log/h10maker_%A.out
#SBATCH --error=log/h10maker_%A.err


export NUM_THREADS=$SLURM_CPUS_PER_TASK
# export NUM_THREADS=2

export INPUT_FOLDER=/cache/clas/e1d/production/pass2/v2/data
export OUTPUT_FOLDER=/work/clas/clase1/tylern/e1d/data/golden_run

mkdir -p $OUTPUT_FOLDER

echo $HOSTNAME

lscpu

echo "STARTING RUN"

STARTTIME=$(date +%s)

cat all.list | parallel -j$NUM_THREADS singularity exec -B /cache:/cache -B /volatile:/volatile -B /work:/work /cvmfs/singularity.opensciencegrid.org/tylern4/clas6:latest h10maker -r $INPUT_FOLDER/{} $OUTPUT_FOLDER/h10_{.}.root

ENDTIME=$(date +%s)
echo "Hostname: $HOSTNAME"
echo "Total runtime: $(($ENDTIME-$STARTTIME))"

echo "DONE!!"
