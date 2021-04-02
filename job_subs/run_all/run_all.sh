#!/bin/bash


echo "Cleanup"
rm -rf /work/tylerns/e1d/outputs/*.root
rm -rf /work/tylerns/e1d/outputs/csv/*
rm -rf /work/tylerns/e1d/outputs/plots/*
#rm -rf /work/tylerns/e1f/outputs/*.root
rm -rf /work/tylerns/e1f/outputs/csv/*
rm -rf /work/tylerns/e1f/outputs/plots/*


RES_E1D_DATA=$(sbatch job_e1d.sh);

echo "Running tests";
RES_E1D_MC=$(sbatch csv_maker_mc_e1d_testing.sh);

#echo "Full Run";
#RES_E1D_MC=$(sbatch csv_maker_mc_e1d.sh);


echo ${RES_E1D_DATA}
echo ${RES_E1D_MC}

# sbatch job_e1d_mc.sh

# sbatch --dependency=afterok:${RES_E1D_DATA##* },${RES_E1D_MC##* } plots_e1d.sh
# sbatch --dependency=afterok:${RES_E1D_DATA##* },${RES_E1D_MC##* } yeilds_e1d.sh
# sbatch --dependency=afterok:${RES_E1D_DATA##* },${RES_E1D_MC##* } cross_sections_e1d.sh
sbatch --dependency=afterok:${RES_E1D_DATA##* },${RES_E1D_MC##* } cross_sections_csv.sh
# sbatch --dependency=afterok:${RES_E1D_DATA##* },${RES_E1D_MC##* } cross_sections_highw.sh


