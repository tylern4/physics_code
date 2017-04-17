#!/bin/bash
source /root/.bashrc
declare -x CLAS_CALDB_DBNAME="calib"
declare -x CLAS_CALDB_HOST=$CLASDB_PORT_3306_TCP_ADDR
declare -x CLAS_CALDB_PASS=""
declare -x CLAS_CALDB_RUNINDEX="RunIndex"
declare -x RECSIS_RUNTIME="/clas/parms/recsis/runtime"
#declare -x CLAS_CALDB_USER=$CLASDB_ENV_MYSQL_USER

function noisyexe {
    # Prints the command to stdout before running
    echo $RED$@$DEF
    $(echo $@)
}

#BASE=/root/code
#DATA=/root/data
DATA=.
BASE=.
VOL=/root

BLUE="\033[34m"
DEF="\033[0m"
MAG="\033[35m"
RED="\033[31m"
GREEN="\033[32m"

echo -e "$RED#################$DEF"
aao_rad < $BASE/aao_rad.inp
part2mctk aao_rad.evt $DATA/thrown.bos
nt10maker -t2 -o$DATA/thrown.hbook $DATA/thrown.bos
h2root $DATA/thrown.hbook $DATA/thrown.root
echo -e "$RED#################$DEF"

echo -e "$BLUE##########$DEF"
gsim_bat -nomcdata -ffread $BASE/gsim.inp -mcin $DATA/thrown.bos -kine 1 -bosout $DATA/gsim.bos
cp $DATA/gsim.bos $DATA/gsim_uncooked.bos
cp $DATA/gsim.bos $DATA/uncooked.bos
noisyexe user_ana -t user_ana.tcl
rm anamonhist logfile
rm $DATA/uncooked.bos
rm $DATA/gsim.bos
mv cooked.bos $DATA/gsim.bos
mv cooked_chist.hbook $DATA/gsim_chist.hbook
nt10maker -t3 $DATA/gsim.bos -o$DATA/gsim.hbook
h2root $DATA/gsim.hbook $DATA/gsim.root
h2root $DATA/gsim_chist.hbook $DATA/gsim_chist.root
echo -e "$BLUE##########$DEF"

echo -e "$GREEN##########$DEF"
export CLAS_CALDB_RUNINDEX=calib_user.RunIndexe1_6
cp gsim_uncooked.bos gpp_in.bos
gpp gsim.inp
rm gpp_in.bos
mv gpp_out.bos reconstructed_GPP_OUT
splitbos -runnum 10 -o reconstructed_gpp_split_out.bos reconstructed_GPP_OUT
cp reconstructed_gpp_split_out.bos uncooked.bos
user_ana -t user_ana.tcl
mv cooked_chist.hbook reconstructed_chist.hbook
mv cooked.bos reconstructed.bos
nt10maker -t3 reconstructed.bos -oreconstructed.hbook
h2root reconstructed.hbook reconstructed.root
h2root reconstructed_chist.hbook reconstructed_chist.root
echo -e "$GREEN##########$DEF"
