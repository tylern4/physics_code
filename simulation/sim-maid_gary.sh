#!/bin/bash
# environment
#source /group/clas/builds/environment.sh
source $HOME/.bashrc
# arg parsing
if [ $# -ne 3 ] ; then
    (echo "Usage: sim-maid.sh <aao_rad_input.inp> <smearbos.inp> <output_prefix>"
     echo ""
     echo "Ensure that TCL scripts fed to user_ana_evan name output correctly;"
     echo "see ~/config/user_ana/e1e.tcl for correct names.") 1>&2
    exit 1
fi
# Existing destination files will be overwritten by mv, so no need to
# handle remove manually.
aao_rad_inp=$1
smearbos_inp=$2
outputprefix=$3
# config files
inp_recsis_e1e_tcl=$HOME/config/user_ana/e1e.tcl
gsim_bat_ffread=$HOME/config/gsim_bat/markov.inp

# PROCESSING

# MAID simulation
$HOME/bin/Linux64RHEL5/aao_rad_gary < $aao_rad_inp
# output is:
# aao_rad.evt: BOS file
# aao_rad.out: unused
# 'aao_rad.rz  p':  unused
# aao_rad.sum: unused
rm aao_rad.out 'aao_rad.rz  p' aao_rad.sum

# Fermi smearing
#smearbos config files are under
#$HOME/config/smearbos/
smearbos < $smearbos_inp
# output is:
# smeared_thrown.bos 

# Detector simulation
gsim_bat -nomcdata -ffread $gsim_bat_ffread -mcin smeared_thrown.bos -kine 1 -bosout GSIM.OUT
# output is:
# GSIM.OUT

# Cooking GSIM output
ln -s -T GSIM.OUT uncooked.bos
user_ana_evan -t $inp_recsis_e1e_tcl
# output is:
# gsim_cooked.bos: cooked BOS output from gpp
# gsim_cooked_chist.hbook: hbook of diagnostic histograms, converted later
# anamonhist: unused
# logfile: unused
# uncooked.bos: unused
rm anamonhist logfile
unlink uncooked.bos

mv cooked.bos gsim_cooked.bos
mv cooked_chist.hbook gsim_cooked_chist.hbook

# Radiative corrections
export CLAS_CALDB_RUNINDEX=calib_user.RunIndexe1_6
gpp  -oGPP_OUT -a2.35 -b2.35 -c2.35 -f0.97 -P0x1b -R36557 ./GSIM.OUT
# output is:
# GPP_OUT: BOS file
# gpp.hbook: unused
rm gpp.hbook
# remove unnecessary GSIM.OUT
rm GSIM.OUT

# I'm not sure what splitbos is doing since -runnum is not a
# documented option.
splitbos -runnum 10 -o gpp_split_out.bos GPP_OUT
# output is:
# gpp_split_out.bos
rm GPP_OUT

# Cleanup
unsetenv CLAS_CALDB_RUNINDEX

# Cook GPP output
ln -s -T gpp_split_out.bos uncooked.bos
user_ana_evan -t $inp_recsis_e1e_tcl
# output is:
# gpp_cooked.bos: cooked BOS output from gpp
# gpp_cooked_chist.hbook: hbook of diagnostic histograms, converted later
# anamonhist: unused
# logfile: unused
# uncooked.bos: unused
rm anamonhist logfile
unlink uncooked.bos
rm gpp_split_out.bos

mv cooked_chist.hbook gpp_cooked_chist.hbook
mv cooked.bos gpp_cooked.bos

# Create ROOT files for thrown and smeared thrown
nt10maker -t2 -oaao_rad.hbook aao_rad.evt
nt10maker_mctk_new -t2 -osmeared_thrown.hbook smeared_thrown.bos
h2root aao_rad.hbook thrown.root
h2root smeared_thrown.hbook smeared_thrown.root
# output is:
# thrown.root: ROOT file with thrown data
# aao_rad.hbook: unused
rm aao_rad.hbook aao_rad.evt
rm smeared_thrown.hbook smeared_thrown.bos

# Create ROOT file for gsim reconstructed
nt10maker_mctk_new -t3 gsim_cooked.bos -ogsim_reconstructed.hbook
h2root gsim_reconstructed.hbook gsim_reconstructed.root
rm gsim_reconstructed.hbook gsim_cooked.bos

# Create ROOT file for reconstructed
nt10maker_mctk_new -t3 gpp_cooked.bos -ogpp_reconstructed.hbook
h2root gpp_reconstructed.hbook gpp_reconstructed.root
# output is:
# reconstructed.root: ROOT file for reconstructed
# reconstructed.hbook: unused
rm gpp_reconstructed.hbook gpp_cooked.bos

# Create ROOT file for chist of gsim reconstructed
h2root gsim_cooked_chist.hbook gsim_reconstructed_chist.root
rm gsim_cooked_chist.hbook

# Create ROOT file for chist of gpp reconstructed
h2root gpp_cooked_chist.hbook gpp_reconstructed_chist.root
rm gpp_cooked_chist.hbook
# output is:
# reconstructed_chist.root

# Final Output are:
# thrown.root
# smeared_thrown.root
# gsim_reconstructed.root
# gsim_reconstructed_chist.root
# gpp_reconstructed.root
# gpp_reconstructed_chist.root

mv thrown.root ${outputprefix}thrown.root
mv smeared_thrown.root ${outputprefix}smeared_thrown.root
mv gsim_reconstructed.root ${outputprefix}gsim_reconstructed.root
mv gsim_reconstructed_chist.root ${outputprefix}gsim_reconstructed_chist.root
mv gpp_reconstructed.root ${outputprefix}gpp_reconstructed.root
mv gpp_reconstructed_chist.root ${outputprefix}gpp_reconstructed_chist.root
