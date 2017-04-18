#!/bin/bash
source /clas-software/env.sh
declare -x CLAS_CALDB_DBNAME="calib"
declare -x CLAS_CALDB_HOST=127.0.0.1
declare -x CLAS_CALDB_PASS=""
declare -x CLAS_CALDB_RUNINDEX="RunIndex"
declare -x RECSIS_RUNTIME="/clas/parms/recsis/runtime"

function noisyexe {
    # Prints the command to stdout before running
    echo $RED$@$DEF
    $(echo $@)
}

CODE=.
#DATA=/srv/data
mkdir -p /srv/data
DATA=/srv/data
BASE=.
VOL=/root

BLUE="\033[34m"
DEF="\033[0m"
MAG="\033[35m"
RED="\033[31m"
GREEN="\033[32m"

echo -e "$RED#################$DEF"
aao_rad < $CODE/aao_rad.inp
mv aao_rad.evt $DATA/aao_rad.evt
part2mctk $DATA/aao_rad.evt $DATA/aao_rad.bos
h10maker -m $DATA/aao_rad.bos $DATA/aao_rad.root

h10maker -p $DATA/aao_rad.evt $DATA/aao_rad_evt_p.root
echo -e "$RED#################$DEF"

echo -e "$BLUE##########$DEF"
gsim_bat -nomcdata -ffread $CODE/gsim.inp -mcin $DATA/aao_rad.bos -kine 1 -bosout $DATA/gsim.bos
h10maker -m $DATA/gsim.bos $DATA/gsim_m.root

gsim_bat -nomcdata -ffread $CODE/gsim.inp -mcin $DATA/aao_rad.evt -kine 1 -bosout $DATA/gsim_evt.bos
h10maker -p $DATA/gsim_evt.bos $DATA/gsim_evt_p.root
#cp gsim.bos uncooked.bos
#noisyexe user_ana -t user_ana.tcl 2> user_ana_1.err 1> user_ana_1.out
#h10maker -r cooked.bos cooked_1.root
#h2root cooked_chist.hbook cooked_chist_1.root
#rm -r uncooked.bos cooked.bos cooked_chist.hbook
#cp gsim_evt.bos uncooked.bos
#noisyexe user_ana -t user_ana.tcl 2> user_ana_2.err 1> user_ana_2.out
#h10maker -r cooked.bos cooked_2.root
#h2root cooked_chist.hbook cooked_chist_2.root
#rm -r uncooked.bos cooked.bos cooked_chist.hbook
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
