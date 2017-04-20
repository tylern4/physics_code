#!/bin/bash
source /clas-software/env.sh
export CLAS_CALDB_DBNAME="calib"
export CLAS_CALDB_PASS=""
export CLAS_CALDB_RUNINDEX="RunIndex"
export RECSIS_RUNTIME="/clas/parms/recsis/runtime"

function noisyexe {
    # Prints the command to stdout before running
    echo $RED$@$DEF
    $(echo $@)
}

BLUE="\033[34m"
DEF="\033[0m"
MAG="\033[35m"
RED="\033[31m"
GREEN="\033[32m"


HOST=${hostname}
if ! [[ $HOST =~ ^[0-9A-F]{6}$ ]]; then
  echo "Using Docker"
  export SYSTEM="Docker"
  export CLAS_CALDB_HOST=$CLASDB_PORT_3306_TCP_ADDR
else
  export SYSTEM="Not Docker"
  export CLAS_CALDB_HOST='127.0.0.1'
fi

CODE=.
mkdir -p data
DATA=./data
BASE=.
VOL=.
mkdir -p logs
LOG=./logs

echo -e "$RED#################$DEF"
aao_rad < $CODE/aao_rad.inp 2> $LOG/aao_rad.err 1> $LOG/aao_rad.out
mv aao_rad.evt $DATA/aao_rad.evt

nt10maker -t2 $DATA/aao_rad.evt -o$DATA/aao_rad.hbook
h2root $DATA/aao_rad.hbook $DATA/aao_rad.root
echo -e "$RED#################$DEF"

echo -e "$BLUE##########$DEF"
gsim_bat -nomcdata -ffread $CODE/gsim.inp -mcin $DATA/aao_rad.evt -kine 1 -bosout $DATA/gsim.evt 2> $LOG/gsim_bat.err 1> $LOG/gsim_bat.out

cp $DATA/gsim.evt uncooked.bos
user_ana -t user_ana.tcl 2> $LOG/user_ana.err 1> $LOG/user_ana.out

cp cooked.bos $DATA/cooked.bos
nt10maker -t1 $DATA/cooked.bos -o$DATA/reconstructed.hbook
h2root $DATA/reconstructed.hbook $DATA/reconstructed.root
echo -e "$BLUE##########$DEF"

echo -e "$GREEN##########$DEF"
##export CLAS_CALDB_RUNINDEX=calib_user.RunIndexe1_6
##cat gpp.inp | grep -v '^#.*' | xargs gpp
#cp gsim_uncooked.bos gpp_in.bos
#gpp gsim.inp
#rm gpp_in.bos
#mv gpp_out.bos reconstructed_GPP_OUT
#splitbos -runnum 10 -o reconstructed_gpp_split_out.bos reconstructed_GPP_OUT
#cp reconstructed_gpp_split_out.bos uncooked.bos
#user_ana -t user_ana.tcl
#mv cooked_chist.hbook reconstructed_chist.hbook
#mv cooked.bos reconstructed.bos
#nt10maker -t3 reconstructed.bos -oreconstructed.hbook
#h2root reconstructed.hbook reconstructed.root
#h2root reconstructed_chist.hbook reconstructed_chist.root
echo -e "$GREEN##########$DEF"
