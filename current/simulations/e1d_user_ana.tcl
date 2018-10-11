source /group/clas/builds/centos7/trunk/reconstruction/recsis/recsis_proc.tcl;
# define packages
turnoff ALL;
global_section off;
turnon seb trk cc tof egn lac user pid;
#
inputfile uncooked.bos;
setc chist_filename histfile;
setc log_file_name logfile;
#
#
set torus_current      3375;
set mini_torus_current 6000;
setc outbanknames(1) "TRGSHEADHEVTEVNTECPBSCPBDCPBCCPBLCPBPARTMCEVTRGSTBID";
#setc outbanknames(1) "TRGSHEADEC  SC  DC0 CC  ECPIECHBPARTTBIDHEVTEVNTDCPBCCPBSCPBECPBCALLTBERTGBITRKSTBTRSCR SCRCECPCCL01LCPBBMPRTDPL";
setc prlink_file_name "prlink_NEW_60_75.bos";
setc bfield_file_name "bgrid_T67to33.fpk";
outputfile cooked.bos PROC 2047;
#
set lseb_nt_do  -1;
set lall_nt_do  -1;
set lscr_nt_do  -1;
#set lmctk_nt_do -1;
set lseb_hist   -1;
set lseb_h_do   -1;
set lmon_hist   -1;
set ltrk_h_do   -1;
set legn_h_do   -1;
set ltof_h_do   -1;
set lec1_h_do   -1;
set lfec_hist   -1;
set lfec_h_do   -1;
set lpart_nt_do -1;
#set lmysql      -1;
#set nmysql      -1;
#
#
fpack "timestop -9999999999"
set lscat $false;
set ldisplay_all $false;
setc rec_prompt "CLASCHEF_recsis> ";
#
setc rec_prompt "CLASCHEF_recsis> ";
go 10000000;
exit_pend;
