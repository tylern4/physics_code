/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef BRANCHES_H
#define BRANCHES_H
#include "TChain.h"
#include "constants.hpp"

/* My Branches */
extern Float_t W;
extern Float_t Q2;
extern Float_t MM;

extern std::vector<bool> *is_proton;
extern std::vector<bool> *is_pip;
extern std::vector<bool> *is_electron;
extern std::vector<bool> *is_pim;

extern std::vector<double> *dt_proton;
extern std::vector<double> *dt_pip;
extern Int_t num_of_pis;
extern bool has_neutron;
//////

// Gary's h10maker
// Declaration of leaf types
extern Int_t npart;
extern Int_t evstat;
extern Int_t intt;
extern Int_t evntid;
extern Int_t evtype;
extern Int_t evntclas;
extern Int_t evthel;
extern Int_t evntclas2;
extern Float_t q_l;
extern Float_t t_l;
extern Float_t tr_time;
extern Float_t rf_time1;
extern Float_t rf_time2;
extern Int_t gpart;
extern Int_t id[MAX_PARTS];   //[gpart]
extern Int_t stat[MAX_PARTS]; //[gpart]
extern Int_t dc[MAX_PARTS];   //[gpart]
extern Int_t cc[MAX_PARTS];   //[gpart]
extern Int_t sc[MAX_PARTS];   //[gpart]
extern Int_t ec[MAX_PARTS];   //[gpart]
extern Int_t lec[MAX_PARTS];  //[gpart]
extern Int_t ccst[MAX_PARTS]; //[gpart]
extern Float_t p[MAX_PARTS];  //[gpart]
extern Float_t m[MAX_PARTS];  //[gpart]
extern Int_t q[MAX_PARTS];    //[gpart]
extern Float_t b[MAX_PARTS];  //[gpart]
extern Float_t cx[MAX_PARTS]; //[gpart]
extern Float_t cy[MAX_PARTS]; //[gpart]
extern Float_t cz[MAX_PARTS]; //[gpart]
extern Float_t vx[MAX_PARTS]; //[gpart]
extern Float_t vy[MAX_PARTS]; //[gpart]
extern Float_t vz[MAX_PARTS]; //[gpart]
extern Int_t dc_part;
extern Int_t dc_sect[MAX_PARTS];   //[dc_part]
extern Int_t dc_trk[MAX_PARTS];    //[dc_part]
extern Int_t dc_stat[MAX_PARTS];   //[dc_part]
extern Float_t dc_vx[MAX_PARTS];   //[dc_part]
extern Float_t dc_vy[MAX_PARTS];   //[dc_part]
extern Float_t dc_vz[MAX_PARTS];   //[dc_part]
extern Float_t dc_vr[MAX_PARTS];   //[dc_part]
extern Float_t dc_xsc[MAX_PARTS];  //[dc_part]
extern Float_t dc_ysc[MAX_PARTS];  //[dc_part]
extern Float_t dc_zsc[MAX_PARTS];  //[dc_part]
extern Float_t dc_cxsc[MAX_PARTS]; //[dc_part]
extern Float_t dc_cysc[MAX_PARTS]; //[dc_part]
extern Float_t dc_czsc[MAX_PARTS]; //[dc_part]
extern Float_t dc_c2[MAX_PARTS];   //[dc_part]
extern Int_t ec_part;
extern Int_t ec_stat[MAX_PARTS]; //[ec_part]
extern Int_t ec_sect[MAX_PARTS]; //[ec_part]
extern Int_t ec_whol[MAX_PARTS]; //[ec_part]
extern Int_t ec_inst[MAX_PARTS]; //[ec_part]
extern Int_t ec_oust[MAX_PARTS]; //[ec_part]
extern Float_t etot[MAX_PARTS];  //[ec_part]
extern Float_t ec_ei[MAX_PARTS]; //[ec_part]
extern Float_t ec_eo[MAX_PARTS]; //[ec_part]
extern Float_t ec_t[MAX_PARTS];  //[ec_part]
extern Float_t ec_r[MAX_PARTS];  //[ec_part]
extern Float_t ech_x[MAX_PARTS]; //[ec_part]
extern Float_t ech_y[MAX_PARTS]; //[ec_part]
extern Float_t ech_z[MAX_PARTS]; //[ec_part]
extern Float_t ec_m2[MAX_PARTS]; //[ec_part]
extern Float_t ec_m3[MAX_PARTS]; //[ec_part]
extern Float_t ec_m4[MAX_PARTS]; //[ec_part]
extern Float_t ec_c2[MAX_PARTS]; //[ec_part]
extern Int_t sc_part;
extern Int_t sc_sect[MAX_PARTS]; //[sc_part]
extern Int_t sc_hit[MAX_PARTS];  //[sc_part]
extern Int_t sc_pd[MAX_PARTS];   //[sc_part]
extern Int_t sc_stat[MAX_PARTS]; //[sc_part]
extern Float_t edep[MAX_PARTS];  //[sc_part]
extern Float_t sc_t[MAX_PARTS];  //[sc_part]
extern Float_t sc_r[MAX_PARTS];  //[sc_part]
extern Float_t sc_c2[MAX_PARTS]; //[sc_part]
extern Int_t cc_part;
extern Int_t cc_sect[MAX_PARTS]; //[cc_part]
extern Int_t cc_hit[MAX_PARTS];  //[cc_part]
extern Int_t cc_segm[MAX_PARTS]; //[cc_part]
extern Int_t nphe[MAX_PARTS];    //[cc_part]
extern Float_t cc_t[MAX_PARTS];  //[cc_part]
extern Float_t cc_r[MAX_PARTS];  //[cc_part]
extern Float_t cc_c2[MAX_PARTS]; //[cc_part]
extern Int_t lac_part;
extern Int_t lec_sect[MAX_PARTS];   //[lac_part]
extern Int_t lec_hit[MAX_PARTS];    //[lac_part]
extern Int_t lec_stat[MAX_PARTS];   //[lac_part]
extern Float_t lec_etot[MAX_PARTS]; //[lac_part]
extern Float_t lec_ein[MAX_PARTS];  //[lac_part]
extern Float_t lec_t[MAX_PARTS];    //[lac_part]
extern Float_t lec_r[MAX_PARTS];    //[lac_part]
extern Float_t lec_x[MAX_PARTS];    //[lac_part]
extern Float_t lec_y[MAX_PARTS];    //[lac_part]
extern Float_t lec_z[MAX_PARTS];    //[lac_part]
extern Float_t lec_c2[MAX_PARTS];   //[lac_part]
extern Int_t st_part;
extern Int_t st_status[MAX_PARTS]; //[st_part]
extern Float_t st_time[MAX_PARTS]; //[st_part]
extern Float_t st_rtrk[MAX_PARTS]; //[st_part]

// List of branches
extern TBranch *b_npart;     //!
extern TBranch *b_evstat;    //!
extern TBranch *b_intt;      //!
extern TBranch *b_evntid;    //!
extern TBranch *b_evtype;    //!
extern TBranch *b_evntclas;  //!
extern TBranch *b_evthel;    //!
extern TBranch *b_evntclas2; //!
extern TBranch *b_q_l;       //!
extern TBranch *b_t_l;       //!
extern TBranch *b_tr_time;   //!
extern TBranch *b_rf_time1;  //!
extern TBranch *b_rf_time2;  //!
extern TBranch *b_gpart;     //!
extern TBranch *b_id;        //!
extern TBranch *b_stat;      //!
extern TBranch *b_dc;        //!
extern TBranch *b_cc;        //!
extern TBranch *b_sc;        //!
extern TBranch *b_ec;        //!
extern TBranch *b_lec;       //!
extern TBranch *b_ccst;      //!
extern TBranch *b_p;         //!
extern TBranch *b_m;         //!
extern TBranch *b_q;         //!
extern TBranch *b_b;         //!
extern TBranch *b_cx;        //!
extern TBranch *b_cy;        //!
extern TBranch *b_cz;        //!
extern TBranch *b_vx;        //!
extern TBranch *b_vy;        //!
extern TBranch *b_vz;        //!
extern TBranch *b_dc_part;   //!
extern TBranch *b_dc_sect;   //!
extern TBranch *b_dc_trk;    //!
extern TBranch *b_dc_stat;   //!
extern TBranch *b_dc_vx;     //!
extern TBranch *b_dc_vy;     //!
extern TBranch *b_dc_vz;     //!
extern TBranch *b_dc_vr;     //!
extern TBranch *b_dc_xsc;    //!
extern TBranch *b_dc_ysc;    //!
extern TBranch *b_dc_zsc;    //!
extern TBranch *b_dc_cxsc;   //!
extern TBranch *b_dc_cysc;   //!
extern TBranch *b_dc_czsc;   //!
extern TBranch *b_dc_c2;     //!
extern TBranch *b_ec_part;   //!
extern TBranch *b_ec_stat;   //!
extern TBranch *b_ec_sect;   //!
extern TBranch *b_ec_whol;   //!
extern TBranch *b_ec_inst;   //!
extern TBranch *b_ec_oust;   //!
extern TBranch *b_etot;      //!
extern TBranch *b_ec_ei;     //!
extern TBranch *b_ec_eo;     //!
extern TBranch *b_ec_t;      //!
extern TBranch *b_ec_r;      //!
extern TBranch *b_ech_x;     //!
extern TBranch *b_ech_y;     //!
extern TBranch *b_ech_z;     //!
extern TBranch *b_ec_m2;     //!
extern TBranch *b_ec_m3;     //!
extern TBranch *b_ec_m4;     //!
extern TBranch *b_ec_c2;     //!
extern TBranch *b_sc_part;   //!
extern TBranch *b_sc_sect;   //!
extern TBranch *b_sc_hit;    //!
extern TBranch *b_sc_pd;     //!
extern TBranch *b_sc_stat;   //!
extern TBranch *b_edep;      //!
extern TBranch *b_sc_t;      //!
extern TBranch *b_sc_r;      //!
extern TBranch *b_sc_c2;     //!
extern TBranch *b_cc_part;   //!
extern TBranch *b_cc_sect;   //!
extern TBranch *b_cc_hit;    //!
extern TBranch *b_cc_segm;   //!
extern TBranch *b_nphe;      //!
extern TBranch *b_cc_t;      //!
extern TBranch *b_cc_r;      //!
extern TBranch *b_cc_c2;     //!
extern TBranch *b_lac_part;  //!
extern TBranch *b_lec_sect;  //!
extern TBranch *b_lec_hit;   //!
extern TBranch *b_lec_stat;  //!
extern TBranch *b_lec_etot;  //!
extern TBranch *b_lec_ein;   //!
extern TBranch *b_lec_t;     //!
extern TBranch *b_lec_r;     //!
extern TBranch *b_lec_x;     //!
extern TBranch *b_lec_y;     //!
extern TBranch *b_lec_z;     //!
extern TBranch *b_lec_c2;    //!
extern TBranch *b_st_part;   //!
extern TBranch *b_st_status; //!
extern TBranch *b_st_time;   //!
extern TBranch *b_st_rtrk;   //!

void getMorebranchs(TTree *myTree);
void getBranches(TTree *myTree);

#endif
