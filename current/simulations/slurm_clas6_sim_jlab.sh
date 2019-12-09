#!/bin/bash

#INPUT Sequences file

#SBATCH -J e1d_sim
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --output=logfiles/sim_%A_%a.out
#SBATCH --error=logfiles/sim_%A_%a.err
#SBATCH --array=1-1000

########=========== Setup Folder and Copy in Input Files ===========########
mkdir -p $HOME/.recsis
touch $HOME/.recsis/recseq.ini
mkdir -p /scratch/${USER}/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}
cd /scratch/${USER}/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}

cp /home/tylern/physics_code/current/simulations/aao_rad.inp aao_rad.inp
cp /home/tylern/physics_code/current/simulations/gsim.inp .
cp /home/tylern/physics_code/current/simulations/user_ana.tcl .

res1=$(date +%s.%N)
########=========== Run Generator ===========########
echo "============ aao_rad ============"
singularity exec \
-B /u/group:/group \
-B /lustre:/lustre \
-B /w/work:/work \
-B /lustre/expphy/volatile:/volatile \
-B /u/home:/home \
--pwd /scratch/${USER}/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID} \
/work/clas/clase1/tylern/clas6.img aao_rad < aao_rad.inp
echo "============ aao_rad ============"

########=========== Run gsim ===========########
echo "============ gsim_bat ============"
singularity exec \
-B /u/group:/group \
-B /lustre:/lustre \
-B /w/work:/work \
-B /lustre/expphy/volatile:/volatile \
-B /u/home:/home \
--pwd /scratch/${USER}/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID} \
/work/clas/clase1/tylern/clas6.img gsim_bat -nomcdata -ffread gsim.inp -mcin aao_rad.evt -bosout gsim.bos
cp gsim.bos gsim_no_gpp.bos
echo "============ gsim_bat ============"

########=========== Run gpp ===========########
echo "============ gpp ============"
singularity exec \
-B /u/group:/group \
-B /lustre:/lustre \
-B /w/work:/work \
-B /lustre/expphy/volatile:/volatile \
-B /u/home:/home \
--pwd /scratch/${USER}/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID} \
/work/clas/clase1/tylern/clas6.img gpp -ouncooked.bos -a2.35 -b2.35 -c2.35 -f0.97 -P0x1b -R23500 gsim.bos
echo "============ gpp ============"

########=========== Run user_ana ===========########
echo "============ user_ana ============"
singularity exec \
-B /u/group:/group \
-B /lustre:/lustre \
-B /w/work:/work \
-B /lustre/expphy/volatile:/volatile \
-B /u/home:/home \
-B $HOME/.recsis:/recsis \
--pwd /scratch/${USER}/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID} \
/work/clas/clase1/tylern/clas6.img user_ana -t user_ana.tcl
echo "============ user_ana ============"


########=========== Run h10maker ===========########
echo "============ h10maker ============"
singularity exec \
-B /u/group:/group \
-B /lustre:/lustre \
-B /w/work:/work \
-B /lustre/expphy/volatile:/volatile \
-B /u/home:/home \
-B $HOME/.recsis:/recsis \
--pwd /scratch/${USER}/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID} \
/work/clas/clase1/tylern/clas6.img h10maker -rpm cooked.bos all.root
echo "============ h10maker ============"

########=========== Copy all the files to Work for output ===========########
xz -9 -v -T4 /scratch/${USER}/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}/gsim_no_gpp.bos
cp /scratch/${USER}/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}/gsim_no_gpp.bos.xz /work/clas/clase1/${USER}/simulations/gsim/gsim_no_gpp_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.bos.xz
cp /scratch/${USER}/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}/all.root /work/clas/clase1/${USER}/simulations/root/e1d_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.root

rm -rf /scratch/${USER}/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}

res2=$(date +%s.%N)
dt=$(echo "$res2 - $res1" | bc)
dd=$(echo "$dt/86400" | bc)
dt2=$(echo "$dt-86400*$dd" | bc)
dh=$(echo "$dt2/3600" | bc)
dt3=$(echo "$dt2-3600*$dh" | bc)
dm=$(echo "$dt3/60" | bc)
ds=$(echo "$dt3-60*$dm" | bc)
echo "Hostname: $HOSTNAME"
printf "Total runtime: %d:%02d:%02d:%02.4f\n" $dd $dh $dm $ds