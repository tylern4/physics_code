#!/bin/bash

#INPUT Sequences file

#SBATCH -J e1d_plotter
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -p BigMem,msmoms
##SBATCH -p BigMem,defq-48core,gpu-v100-16gb,gpu-v100-32gb
#SBATCH --output=log/e1d_plotter_%A.out
#SBATCH --error=log/e1d_plotter_%A.err

module load gcc/6.4.0
#module load python3/anaconda/5.2.0
module load python3/anaconda/2020.02

export CC=$(which gcc)
export CXX=$(which g++)


#CERN ROOT
export ROOTSYS=/work/gothelab/software/root
export PATH=$ROOTSYS/bin:$PATH
export PYTHONDIR=$ROOTSYS
export LD_LIBRARY_PATH=$ROOTSYS/lib:$PYTHONDIR/lib:$ROOTSYS/bindings/pyroot:$LD_LIBRARY_PATH
export PYTHONPATH=/usr/local/lib:$ROOTSYS/lib:$PYTHONPATH:$ROOTSYS/bindings/pyroot

export CLAS_PARMS=/work/gothelab/clas6/parms

export MPLBACKEND="module://gr.matplotlib.backend_gr"

cd /home/tylerns/physics_code/python

echo "plots_e1d"
rm -rf /work/tylerns/e1d/outputs/plots
mkdir -p /work/tylerns/e1d/outputs/plots
/work/apps/python3/anaconda3/2020.02/bin/python acceptance.py --mc /work/tylerns/e1d/outputs/csv/mc_osg_e1d.csv --data /work/tylerns/e1d/outputs/csv/data_e1d.feather --folder /work/tylerns/e1d/outputs/plots


