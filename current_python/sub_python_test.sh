#!/bin/sh
#SBATCH --job-name=e1d_python
#SBATCH --output /home/tylerns/job_sub/logfiles/e1dpy_%j.out
#SBATCH --error /home/tylerns/job_sub/logfiles/e1dpy_%j.err
#SBATCH -p all
###Number of Cores Max 20
#SBATCH -n 1
###Number of Nodes
#SBATCH -N 1

source /share/apps/Modules/3.2.10/init/modules.sh
module load root/6.04-minuit
module load gcc/5.4.0
module load python/anaconda2.3

export PYTHONPATH=$PYTHONPATH:$(root-config --libdir)
export PATH=$PATH:/home/tylerns/.local/bin


###script command here
python /home/tylerns/physics_code/current_python/test.py