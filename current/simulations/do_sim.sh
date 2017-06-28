#!/bin/bash
source /clas-software/env.sh
export CLAS_CALDB_DBNAME="calib"
export CLAS_CALDB_PASS=""
export CLAS_CALDB_RUNINDEX="RunIndex"
export RECSIS_RUNTIME="/clas/parms/recsis/runtime"

export SYSTEM="Docker"
export CLAS_CALDB_HOST=$CLASDB_PORT_3306_TCP_ADDR

aao_rad < aao_rad.inp >/dev/null 2>&1

gsim_bat -nomcdata -ffread gsim.inp -mcin aao_rad.evt -kine 1 -bosout uncooked.bos >/dev/null 2>&1

user_ana -t user_ana.tcl >/dev/null 2>&1
nt10maker -t1 cooked.bos -ocooked.hbook >/dev/null 2>&1
h2root cooked.hbook cooked.root >/dev/null 2>&1
