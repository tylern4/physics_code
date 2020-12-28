#!/bin/bash

#INPUT Sequences file

#SBATCH -J e1f_plotter
#SBATCH -n 8
#SBATCH -N 1
#SBATCH -p msmoms,BigMem
#SBATCH --output=log/e1f_plotter_%A.out
#SBATCH --error=log/e1f_plotter_%A.err

module load gcc/6.4.0
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

echo "plots_e1f"
rm -rf /work/tylerns/e1f/outputs/csv/plots
mkdir -p /work/tylerns/e1f/outputs/csv/plots
/work/apps/python3/anaconda3/2020.02/bin/python acceptance.py --e1f --mc /work/tylerns/e1f/outputs/csv/mc_osg_e1f.csv --data /work/tylerns/e1f/outputs/csv/data_e1f.feather --folder /work/tylerns/e1f/outputs/plots

