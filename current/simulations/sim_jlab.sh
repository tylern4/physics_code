#!/bin/bash
source /site/12gev_phys/softenv.sh 2.2
export CLAS6=/work/clas/clase1/tylern/clas-software/build
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CLAS6/lib
export PATH=$PATH:$CLAS6/bin
export CLAS_PARMS=/group/clas/parms

aao_rad < aao_rad.inp

gsim_bat -nomcdata -ffread gsim.inp -mcin aao_rad.evt -bosout gsim.bos #>/dev/null 2>&1

gpp -ouncooked.bos -a2.35 -b2.35 -c2.35 -f0.97 -P0x1b -R36557 gsim.bos #>/dev/null 2>&1

user_ana -t user_ana.tcl #>/dev/null 2>&1


h10maker -r cooked.bos cooked.root #>/dev/null 2>&1
h10maker -p cooked.bos cooked_mc.root #>/dev/null 2>&1
gzip gsim.bos
gzip aao_rad.evt
gzip cooked.bos
