#!/bin/bash

RES_E1D_DATA=$(sbatch job_e1d.sh);
RES_E1F_DATA=$(sbatch job_e1f.sh);

RES_MC=$(sbatch csv_maker_mc.sh);
#RES_MC_OLD=$(sbatch csv_maker_mc_old.sh);

echo ${RES_E1D_DATA}
echo ${RES_E1F_DATA}


sbatch --dependency=afterok:${RES_E1F_DATA##* },${RES_E1D_DATA##* },${RES_MC##* } plots.sh



