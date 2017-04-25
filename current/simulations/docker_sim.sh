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

DATE=`date +%m-%d-%y_%H-%M-%S`
HOST=`hostname`


if ! [[ $HOST =~ ^[0-9A-F]{6}$ ]]; then
  echo -e $MAG"Using Docker!!!!!!!!!"$DEF
  echo -e $GREEN"$HOST"$DEF
  export SYSTEM="Docker"
  export CLAS_CALDB_HOST=$CLASDB_PORT_3306_TCP_ADDR
  CODE=/root/code
  DATA=/root/data/out_$DATE
  VOL=/root/
  LOG=$DATA/logs
  mkdir -p $DATA
  mkdir -p $LOG
else
  export SYSTEM="Not Docker"
  export CLAS_CALDB_HOST='127.0.0.1'
  CODE=.
  mkdir -p ./data
  DATA=./data
  VOL=.
  mkdir -p ./data/logs
  LOG=./data/logs
fi

evt2root () {
  TYPE=$1
  FILE=$2
  BNAME=${FILE%.evt}
  nt10maker -t$TYPE $FILE -o$VOL/temp_$RANDOM.hbook 1> $VOL/log.out
  h2root $VOL/temp_$RANDOM.hbook $BNAME.root 1> $VOL/log.out
}

bos2root () {
  TYPE=$1
  FILE=$2
  BNAME=${FILE%.bos}
  nt10maker -t$TYPE $FILE -o$VOL/temp_$RANDOM.hbook 1> $VOL/log.out
  h2root $VOL/temp_$RANDOM.hbook $BNAME.root 1> $VOL/log.out
  echo $BNAME
}

echo -e "$RED####\tAAO_RAD\t####$DEF"
aao_rad < $CODE/aao_rad.inp 2> $LOG/aao_rad.err 1> $LOG/aao_rad.out
mv aao_rad.out $LOG/aao_rad_$DATE.out
mv aao_rad.evt $DATA/aao_rad_$DATE.evt
rm -f *.rz* *.sum

nt10maker -t2 $DATA/aao_rad_$DATE.evt -o$VOL/temp_$RANDOM.hbook 1> $VOL/log.out
h2root $VOL/temp_$RANDOM.hbook aao_rad_$DATE.root 1> $VOL/log.out
echo -e "$RED####\tAAO_RAD\t####$DEF"



echo -e "$BLUE####\tGSIM\t####$DEF"
if [[ $SYSTEM == "Docker" ]]; then
  echo -e $MAG"Using Docker!!!!!!!!!"$DEF
  splitbos -o $DATA/thrown_$DATE -seq 1 -q -n 10000 -M 8 $DATA/aao_rad_$DATE.evt
  parallel --will-cite --jobs 8 "gsim_bat -nomcdata -ffread $CODE/gsim.inp -mcin $DATA/thrown_$DATE.0{#} -kine 1 -bosout $DATA/gsim_$DATE.0{#}.evt 2> $LOG/gsim_bat_0{#}.err 1> $LOG/gsim_bat_0{#}.out" ::: {0..7}
  parallel --will-cite --jobs 1 "cp $DATA/gsim_$DATE.0{#}.evt uncooked.bos && user_ana -t user_ana.tcl 2> $LOG/user_ana_0{#}.err 1> $LOG/user_ana_0{#}.out && mv cooked.bos $DATA/cooked_$DATE.0{#}.bos" ::: {0..7}
  parallel --will-cite --jobs 8 "nt10maker -t1 $DATA/cooked_$DATE.0{#}.bos -o$DATA/cooked_$DATE.0{#}.hbook 2> $LOG/nt10_0{#}.err 1> $LOG/nt10_0{#}.out" ::: {0..7}
  parallel --will-cite --jobs 8 "h2root $DATA/cooked_$DATE.0{#}.hbook $DATA/cooked_$DATE.0{#}.root 2> $LOG/h2root_0{#}.err 1> $LOG/h2root_0{#}.out" ::: {0..7}
else
  gsim_bat -nomcdata -ffread $CODE/gsim.inp -mcin $DATA/aao_rad.evt -kine 1 -bosout $DATA/gsim_$DATE.evt 2> $LOG/gsim_bat.err 1> $LOG/gsim_bat.out
  cp $DATA/gsim.evt uncooked.bos
  cp $DATA/gsim_$DATE.evt uncooked.bos
  user_ana -t user_ana.tcl 2> $LOG/user_ana.err 1> $LOG/user_ana.out
  mv cooked.bos $DATA/cooked_$DATE.bos
  nt10maker -t1 $DATA/cooked_$DATE.bos -o$DATA/cooked_$DATE.hbook
  h2root $DATA/cooked_$DATE.hbook $DATA/cooked_$DATE.root
fi
echo -e "$BLUE####\tGSIM\t####$DEF"


#echo -e "$GREEN####\tGPP\t####$DEF"
#export CLAS_CALDB_RUNINDEX=calib_user.RunIndexe1_6
#if [[ $SYSTEM == "Docker" ]]; then
#  echo -e $MAG"Using Docker!!!!!!!!!"$DEF
#  parallel --will-cite --jobs 1 "cp $DATA/cooked_$DATE.0{#}.bos gpp_in.bos && cat $CODE/gpp.inp | grep -v '^#.*' | xargs gpp 2> $LOG/gpp_0{#}.err 1> $LOG/gpp_0{#}.out && mv gpp_out.bos $DATA/gpp_$DATE.0{#}.bos" ::: {0..7}
#  parallel --will-cite --jobs 1 "cp $DATA/gpp_$DATE.0{#}.bos uncooked.bos && user_ana -t user_ana.tcl 2> $LOG/user_ana_gpp_0{#}.err 1> $LOG/user_ana_gpp_0{#}.out && mv cooked.bos $DATA/gpp_cooked_$DATE.0{#}.bos" ::: {0..7}
#  parallel --will-cite --jobs 8 "nt10maker -t1 $DATA/gpp_cooked_$DATE.0{#}.bos -o$DATA/gpp_cooked_$DATE.0{#}.hbook 2> $LOG/nt10_0{#}.err 1> $LOG/nt10_0{#}.out" ::: {0..7}
#  parallel --will-cite --jobs 8 "h2root $DATA/gpp_cooked_$DATE.0{#}.hbook $DATA/gpp_cooked_$DATE.0{#}.root 2> $LOG/h2root_0{#}.err 1> $LOG/h2root_0{#}.out" ::: {0..7}
#else
#  echo "NOT SETUP YET!"
#fi
#echo -e "$GREEN####\tGPP\t####$DEF"



echo -e "$RED####\tCLEANUP\t####$DEF"
rm -f *.evt aao_rad.out aao_rad.rz* aao_rad.sum *.hbook *.bos *.root *.log *.rzn logfile anamonhist
echo -e "$GREEN####\tDONE\t####$DEF"
