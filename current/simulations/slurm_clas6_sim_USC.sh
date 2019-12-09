#!/bin/bash

#INPUT Sequences file

#SBATCH -J e1d_sim
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --output=logfiles/sim_%A_%a.out
#SBATCH --error=logfiles/sim_%A_%a.err
#SBATCH --array=1-1000

export SING_IMG=/work/gothelab/clas6/clas6_new.img
export DATA_DIR=/work/gothelab/clas6

########=========== Setup Folder and Copy in Input Files ===========########
mkdir -p $HOME/.recsis
touch $HOME/.recsis/recseq.ini
mkdir -p /local/${USER}/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}
cd /local/${USER}/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}

cp $HOME/physics_code/current/simulations/aao_rad_test.inp .
cp $HOME/physics_code/current/simulations/gsim.inp .
cp $HOME/physics_code/current/simulations/user_ana.tcl .

res1=$(date +%s.%N)
########=========== Run Generator ===========########
echo "============ aao_rad ============"
singularity exec \
-B /work:/work \
-B ${DATA_DIR}:/data \
--pwd /local/${USER}/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID} \
${SING_IMG} aao_rad < aao_rad.inp
echo "============ aao_rad ============"

########=========== Run gsim ===========########
echo "============ gsim_bat ============"
singularity exec \
-B /work:/work \
-B ${DATA_DIR}:/data \
--pwd /local/${USER}/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID} \
${SING_IMG} gsim_bat -nomcdata -ffread gsim.inp -mcin aao_rad.evt -bosout gsim.bos
cp gsim.bos gsim_no_gpp.bos
echo "============ gsim_bat ============"

########=========== Run gpp ===========########
echo "============ gpp ============"
singularity exec \
-B /work:/work \
-B ${DATA_DIR}:/data \
--pwd /local/${USER}/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID} \
${SING_IMG} gpp -ouncooked.bos -a2.35 -b2.35 -c2.35 -f0.97 -P0x1b -R23500 gsim.bos
echo "============ gpp ============"

########=========== Run user_ana ===========########
echo "============ user_ana ============"
singularity exec \
-B /work:/work \
-B ${DATA_DIR}:/data \
-B $HOME/.recsis:/recsis \
--pwd /local/${USER}/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID} \
${SING_IMG} user_ana -t user_ana.tcl
echo "============ user_ana ============"


########=========== Run h10maker ===========########
echo "============ h10maker ============"
singularity exec \
-B /work:/work \
-B ${DATA_DIR}:/data \
--pwd /local/${USER}/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID} \
${SING_IMG} h10maker -rpm cooked.bos all.root
echo "============ h10maker ============"

########=========== Copy all the files to Work for output ===========########
xz -9 -v -T4 /local/${USER}/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}/gsim_no_gpp.bos
cp /local/${USER}/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}/gsim_no_gpp.bos.xz ${DATA_DIR}/e1d/sim/gsim/gsim_no_gpp_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.bos.xz
cp /local/${USER}/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}/all.root ${DATA_DIR}/e1d/sim/root/e1d_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.root

rm -rf /local/${USER}/${SLURM_JOB_ID}/${SLURM_ARRAY_TASK_ID}

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