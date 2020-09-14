#!/bin/bash

RES_E1D_DATA=$(sbatch job_e1d.sh);
RES_E1F_DATA=$(sbatch job_e1f.sh);

RES_E1D_MC=$(sbatch csv_maker_mc_e1d.sh);
RES_E1F_MC=$(sbatch csv_maker_mc_e1f.sh);

echo ${RES_E1D_DATA}
echo ${RES_E1F_DATA}
echo ${RES_E1D_MC}
echo ${RES_E1F_MC}

sbatch job_e1f_mc.sh
sbatch job_e1d_mc.sh


sbatch --dependency=afterok:${RES_E1D_DATA##* },${RES_E1D_MC##* } plots_e1d.sh
sbatch --dependency=afterok:${RES_E1F_DATA##* },${RES_E1F_MC##* } plots_e1f.sh



