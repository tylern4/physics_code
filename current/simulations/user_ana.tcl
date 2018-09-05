source /clas/parms/recsis/recsis_proc.tcl;

# define packages
turnoff ALL;
global_section off;
turnon seb trk cc tof egn lac user pid;
#
inputfile uncooked.bos;
setc chist_filename cooked_chist.hbook;
setc log_file_name logfile;

#
#
set torus_current      3375;
set mini_torus_current 5996;
setc outbanknames(1) "HEADHEVTEVNTECPBSCPBDCPBCCPBLCPBPARTMCEVTRGSTBIDMCTKMCVX";
setc prlink_file_name "prlink_NEW_60_75.bos";
setc bfield_file_name "bgrid_T67to33.fpk";
outputfile cooked.bos PROC 2047;
#



# Franz's tcl variables
#/set trk_maxiter       8;
#/set trk_minhits(1)    2;
#/set trk_lrambfit_chi2 50.;
#/set trk_tbtfit_chi2   70.;
#/set trk_prfit_chi2    70.;
#/set trk_statistics    3 ;
#
#nachalo
#/set ltrk_h_do    -1;
#/set legn_h_do    -1;
#/set lcc_h_do     -1;
#/set ltof_h_do    -1;
#/set lst_h_do      0;
#/set ltime_h_do   -1;
#/set lec1_h_do    -1;
#/set lusr0_h_do    0;
#/set lusr1_h_do    0;

#/set lseb_h_do    -1;
#set lmctk_nt_do     -1;
#/set lseb_hist       -1;
#/set lmon_hist         0;
#/set lfec_hist       -1;
#/set lmvrt_nt_do     -1;




set dc_xvst_choice 0;
set ltrk_do    -1;
set legn_do    -1;
set lcc_do     -1;
set ltof_do    -1;
set lst_do      0;
set ltime_do   0;
set lec1_do    -1;
#1
set ltagger_do 0;
set lseb_do    -1;
#2,3
set lusr0_do   -1;
set lusr1_do   -1;
set lrf_do      0;
set lall_nt_do      -1;
#4,5,6,7
set lechb_nt_do 0;
set ltrk_nt_do 0;
set lec_nt_do 0 ;
set llac_nt_do 0;
set lccr_nt_do  0;
set lscr_nt_do  0;
set ltbtk_nt_do 0;
set lseb_nt_do      -1;
set lpart_nt_do     -1;
set nt_id_cut 0;
set def_geom -1;
set def_adc -1;
set def_tdc -1;
set def_atten  0;
set whole_surf 16.;
set inner_surf 1.;
set outer_surf 28.;
set m2_ecel_cut 1000.;
set m3_ecel_cut 10000.;
set trkec__match 50;

#setc outbanknames(1) "HEADHLS DC0 CC  SC  EC1 EC  CALLPARTTBIDHEVTEVNTDCPBCCPBSCPBECPBLCPBTGBIFBPMCL01TBLATRGSTRPBTBERMCTK";

fpack "timestop -9999999999"
set lscat $false;
set ldisplay_all $false;
#/setc rec_prompt "[exec whoami]_recsis> ";

go 10000000;
#/status;
exit_pend;
