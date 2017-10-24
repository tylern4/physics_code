/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#include "branches.hpp"

/* My Branches */
Float_t W;
Float_t Q2;
Float_t MM;

std::vector<bool> *is_proton;
std::vector<bool> *is_pip;
std::vector<bool> *is_electron;
std::vector<bool> *is_pim;

std::vector<double> *dt_proton = 0;
std::vector<double> *dt_pip = 0;
Int_t num_of_pis;
bool has_neutron;
//////
// Gary's h10maker
Int_t npart;
Int_t evstat;
Int_t intt;
Int_t evntid;
Int_t evtype;
Int_t evntclas;
Int_t evthel;
Int_t evntclas2;
Float_t q_l;
Float_t t_l;
Float_t tr_time;
Float_t rf_time1;
Float_t rf_time2;
Int_t gpart;
Int_t id[MAX_PARTS];   //[gpart]
Int_t stat[MAX_PARTS]; //[gpart]
Int_t dc[MAX_PARTS];   //[gpart]
Int_t cc[MAX_PARTS];   //[gpart]
Int_t sc[MAX_PARTS];   //[gpart]
Int_t ec[MAX_PARTS];   //[gpart]
Int_t lec[MAX_PARTS];  //[gpart]
Int_t ccst[MAX_PARTS]; //[gpart]
Float_t p[MAX_PARTS];  //[gpart]
Float_t m[MAX_PARTS];  //[gpart]
Int_t q[MAX_PARTS];    //[gpart]
Float_t b[MAX_PARTS];  //[gpart]
Float_t cx[MAX_PARTS]; //[gpart]
Float_t cy[MAX_PARTS]; //[gpart]
Float_t cz[MAX_PARTS]; //[gpart]
Float_t vx[MAX_PARTS]; //[gpart]
Float_t vy[MAX_PARTS]; //[gpart]
Float_t vz[MAX_PARTS]; //[gpart]
Int_t dc_part;
Int_t dc_sect[MAX_PARTS];   //[dc_part]
Int_t dc_trk[MAX_PARTS];    //[dc_part]
Int_t dc_stat[MAX_PARTS];   //[dc_part]
Float_t dc_vx[MAX_PARTS];   //[dc_part]
Float_t dc_vy[MAX_PARTS];   //[dc_part]
Float_t dc_vz[MAX_PARTS];   //[dc_part]
Float_t dc_vr[MAX_PARTS];   //[dc_part]
Float_t dc_xsc[MAX_PARTS];  //[dc_part]
Float_t dc_ysc[MAX_PARTS];  //[dc_part]
Float_t dc_zsc[MAX_PARTS];  //[dc_part]
Float_t dc_cxsc[MAX_PARTS]; //[dc_part]
Float_t dc_cysc[MAX_PARTS]; //[dc_part]
Float_t dc_czsc[MAX_PARTS]; //[dc_part]
Float_t dc_c2[MAX_PARTS];   //[dc_part]
Int_t ec_part;
Int_t ec_stat[MAX_PARTS]; //[ec_part]
Int_t ec_sect[MAX_PARTS]; //[ec_part]
Int_t ec_whol[MAX_PARTS]; //[ec_part]
Int_t ec_inst[MAX_PARTS]; //[ec_part]
Int_t ec_oust[MAX_PARTS]; //[ec_part]
Float_t etot[MAX_PARTS];  //[ec_part]
Float_t ec_ei[MAX_PARTS]; //[ec_part]
Float_t ec_eo[MAX_PARTS]; //[ec_part]
Float_t ec_t[MAX_PARTS];  //[ec_part]
Float_t ec_r[MAX_PARTS];  //[ec_part]
Float_t ech_x[MAX_PARTS]; //[ec_part]
Float_t ech_y[MAX_PARTS]; //[ec_part]
Float_t ech_z[MAX_PARTS]; //[ec_part]
Float_t ec_m2[MAX_PARTS]; //[ec_part]
Float_t ec_m3[MAX_PARTS]; //[ec_part]
Float_t ec_m4[MAX_PARTS]; //[ec_part]
Float_t ec_c2[MAX_PARTS]; //[ec_part]
Int_t sc_part;
Int_t sc_sect[MAX_PARTS]; //[sc_part]
Int_t sc_hit[MAX_PARTS];  //[sc_part]
Int_t sc_pd[MAX_PARTS];   //[sc_part]
Int_t sc_stat[MAX_PARTS]; //[sc_part]
Float_t edep[MAX_PARTS];  //[sc_part]
Float_t sc_t[MAX_PARTS];  //[sc_part]
Float_t sc_r[MAX_PARTS];  //[sc_part]
Float_t sc_c2[MAX_PARTS]; //[sc_part]
Int_t cc_part;
Int_t cc_sect[MAX_PARTS]; //[cc_part]
Int_t cc_hit[MAX_PARTS];  //[cc_part]
Int_t cc_segm[MAX_PARTS]; //[cc_part]
Int_t nphe[MAX_PARTS];    //[cc_part]
Float_t cc_t[MAX_PARTS];  //[cc_part]
Float_t cc_r[MAX_PARTS];  //[cc_part]
Float_t cc_c2[MAX_PARTS]; //[cc_part]
Int_t lac_part;
Int_t lec_sect[MAX_PARTS];   //[lac_part]
Int_t lec_hit[MAX_PARTS];    //[lac_part]
Int_t lec_stat[MAX_PARTS];   //[lac_part]
Float_t lec_etot[MAX_PARTS]; //[lac_part]
Float_t lec_ein[MAX_PARTS];  //[lac_part]
Float_t lec_t[MAX_PARTS];    //[lac_part]
Float_t lec_r[MAX_PARTS];    //[lac_part]
Float_t lec_x[MAX_PARTS];    //[lac_part]
Float_t lec_y[MAX_PARTS];    //[lac_part]
Float_t lec_z[MAX_PARTS];    //[lac_part]
Float_t lec_c2[MAX_PARTS];   //[lac_part]
Int_t st_part;
Int_t st_status[MAX_PARTS]; //[st_part]
Float_t st_time[MAX_PARTS]; //[st_part]
Float_t st_rtrk[MAX_PARTS]; //[st_part]

// List of branches
TBranch *b_npart;     //!
TBranch *b_evstat;    //!
TBranch *b_intt;      //!
TBranch *b_evntid;    //!
TBranch *b_evtype;    //!
TBranch *b_evntclas;  //!
TBranch *b_evthel;    //!
TBranch *b_evntclas2; //!
TBranch *b_q_l;       //!
TBranch *b_t_l;       //!
TBranch *b_tr_time;   //!
TBranch *b_rf_time1;  //!
TBranch *b_rf_time2;  //!
TBranch *b_gpart;     //!
TBranch *b_id;        //!
TBranch *b_stat;      //!
TBranch *b_dc;        //!
TBranch *b_cc;        //!
TBranch *b_sc;        //!
TBranch *b_ec;        //!
TBranch *b_lec;       //!
TBranch *b_ccst;      //!
TBranch *b_p;         //!
TBranch *b_m;         //!
TBranch *b_q;         //!
TBranch *b_b;         //!
TBranch *b_cx;        //!
TBranch *b_cy;        //!
TBranch *b_cz;        //!
TBranch *b_vx;        //!
TBranch *b_vy;        //!
TBranch *b_vz;        //!
TBranch *b_dc_part;   //!
TBranch *b_dc_sect;   //!
TBranch *b_dc_trk;    //!
TBranch *b_dc_stat;   //!
TBranch *b_dc_vx;     //!
TBranch *b_dc_vy;     //!
TBranch *b_dc_vz;     //!
TBranch *b_dc_vr;     //!
TBranch *b_dc_xsc;    //!
TBranch *b_dc_ysc;    //!
TBranch *b_dc_zsc;    //!
TBranch *b_dc_cxsc;   //!
TBranch *b_dc_cysc;   //!
TBranch *b_dc_czsc;   //!
TBranch *b_dc_c2;     //!
TBranch *b_ec_part;   //!
TBranch *b_ec_stat;   //!
TBranch *b_ec_sect;   //!
TBranch *b_ec_whol;   //!
TBranch *b_ec_inst;   //!
TBranch *b_ec_oust;   //!
TBranch *b_etot;      //!
TBranch *b_ec_ei;     //!
TBranch *b_ec_eo;     //!
TBranch *b_ec_t;      //!
TBranch *b_ec_r;      //!
TBranch *b_ech_x;     //!
TBranch *b_ech_y;     //!
TBranch *b_ech_z;     //!
TBranch *b_ec_m2;     //!
TBranch *b_ec_m3;     //!
TBranch *b_ec_m4;     //!
TBranch *b_ec_c2;     //!
TBranch *b_sc_part;   //!
TBranch *b_sc_sect;   //!
TBranch *b_sc_hit;    //!
TBranch *b_sc_pd;     //!
TBranch *b_sc_stat;   //!
TBranch *b_edep;      //!
TBranch *b_sc_t;      //!
TBranch *b_sc_r;      //!
TBranch *b_sc_c2;     //!
TBranch *b_cc_part;   //!
TBranch *b_cc_sect;   //!
TBranch *b_cc_hit;    //!
TBranch *b_cc_segm;   //!
TBranch *b_nphe;      //!
TBranch *b_cc_t;      //!
TBranch *b_cc_r;      //!
TBranch *b_cc_c2;     //!
TBranch *b_lac_part;  //!
TBranch *b_lec_sect;  //!
TBranch *b_lec_hit;   //!
TBranch *b_lec_stat;  //!
TBranch *b_lec_etot;  //!
TBranch *b_lec_ein;   //!
TBranch *b_lec_t;     //!
TBranch *b_lec_r;     //!
TBranch *b_lec_x;     //!
TBranch *b_lec_y;     //!
TBranch *b_lec_z;     //!
TBranch *b_lec_c2;    //!
TBranch *b_st_part;   //!
TBranch *b_st_status; //!
TBranch *b_st_time;   //!
TBranch *b_st_rtrk;   //!

void getMorebranchs(TTree *myTree) {
  myTree->SetBranchAddress("W", &W);
  myTree->SetBranchAddress("Q2", &Q2);
  myTree->SetBranchAddress("MM", &MM);

  myTree->SetBranchAddress("is_electron", &is_electron);
  myTree->SetBranchAddress("is_proton", &is_proton);
  myTree->SetBranchAddress("is_pip", &is_pip);
  myTree->SetBranchAddress("is_pim", &is_pim);

  myTree->SetBranchAddress("DeltaT_P", &dt_proton);
  myTree->SetBranchAddress("DeltaT_Pip", &dt_pip);
  myTree->SetBranchAddress("NumPI", &num_of_pis);
  myTree->SetBranchAddress("m", &m);
  myTree->SetBranchStatus("*", 1);
}

void getBranches(TTree *myTree) {
  myTree->SetBranchAddress("npart", &npart, &b_npart);
  myTree->SetBranchAddress("evstat", &evstat, &b_evstat);
  myTree->SetBranchAddress("intt", &intt, &b_intt);
  myTree->SetBranchAddress("evntid", &evntid, &b_evntid);
  myTree->SetBranchAddress("evtype", &evtype, &b_evtype);
  myTree->SetBranchAddress("evntclas", &evntclas, &b_evntclas);
  myTree->SetBranchAddress("evthel", &evthel, &b_evthel);
  myTree->SetBranchAddress("evntclas2", &evntclas2, &b_evntclas2);
  myTree->SetBranchAddress("q_l", &q_l, &b_q_l);
  myTree->SetBranchAddress("t_l", &t_l, &b_t_l);
  myTree->SetBranchAddress("tr_time", &tr_time, &b_tr_time);
  myTree->SetBranchAddress("rf_time1", &rf_time1, &b_rf_time1);
  myTree->SetBranchAddress("rf_time2", &rf_time2, &b_rf_time2);
  myTree->SetBranchAddress("gpart", &gpart, &b_gpart);
  myTree->SetBranchAddress("id", id, &b_id);
  myTree->SetBranchAddress("stat", stat, &b_stat);
  myTree->SetBranchAddress("dc", dc, &b_dc);
  myTree->SetBranchAddress("cc", cc, &b_cc);
  myTree->SetBranchAddress("sc", sc, &b_sc);
  myTree->SetBranchAddress("ec", ec, &b_ec);
  myTree->SetBranchAddress("lec", lec, &b_lec);
  myTree->SetBranchAddress("ccst", ccst, &b_ccst);
  myTree->SetBranchAddress("p", p, &b_p);
  myTree->SetBranchAddress("m", m, &b_m);
  myTree->SetBranchAddress("q", q, &b_q);
  myTree->SetBranchAddress("b", b, &b_b);
  myTree->SetBranchAddress("cx", cx, &b_cx);
  myTree->SetBranchAddress("cy", cy, &b_cy);
  myTree->SetBranchAddress("cz", cz, &b_cz);
  myTree->SetBranchAddress("vx", vx, &b_vx);
  myTree->SetBranchAddress("vy", vy, &b_vy);
  myTree->SetBranchAddress("vz", vz, &b_vz);
  myTree->SetBranchAddress("dc_part", &dc_part, &b_dc_part);
  myTree->SetBranchAddress("dc_sect", dc_sect, &b_dc_sect);
  myTree->SetBranchAddress("dc_trk", dc_trk, &b_dc_trk);
  myTree->SetBranchAddress("dc_stat", dc_stat, &b_dc_stat);
  myTree->SetBranchAddress("dc_vx", dc_vx, &b_dc_vx);
  myTree->SetBranchAddress("dc_vy", dc_vy, &b_dc_vy);
  myTree->SetBranchAddress("dc_vz", dc_vz, &b_dc_vz);
  myTree->SetBranchAddress("dc_vr", dc_vr, &b_dc_vr);
  myTree->SetBranchAddress("dc_xsc", dc_xsc, &b_dc_xsc);
  myTree->SetBranchAddress("dc_ysc", dc_ysc, &b_dc_ysc);
  myTree->SetBranchAddress("dc_zsc", dc_zsc, &b_dc_zsc);
  myTree->SetBranchAddress("dc_cxsc", dc_cxsc, &b_dc_cxsc);
  myTree->SetBranchAddress("dc_cysc", dc_cysc, &b_dc_cysc);
  myTree->SetBranchAddress("dc_czsc", dc_czsc, &b_dc_czsc);
  myTree->SetBranchAddress("dc_c2", dc_c2, &b_dc_c2);
  myTree->SetBranchAddress("ec_part", &ec_part, &b_ec_part);
  myTree->SetBranchAddress("ec_stat", ec_stat, &b_ec_stat);
  myTree->SetBranchAddress("ec_sect", ec_sect, &b_ec_sect);
  myTree->SetBranchAddress("ec_whol", ec_whol, &b_ec_whol);
  myTree->SetBranchAddress("ec_inst", ec_inst, &b_ec_inst);
  myTree->SetBranchAddress("ec_oust", ec_oust, &b_ec_oust);
  myTree->SetBranchAddress("etot", etot, &b_etot);
  myTree->SetBranchAddress("ec_ei", ec_ei, &b_ec_ei);
  myTree->SetBranchAddress("ec_eo", ec_eo, &b_ec_eo);
  myTree->SetBranchAddress("ec_t", ec_t, &b_ec_t);
  myTree->SetBranchAddress("ec_r", ec_r, &b_ec_r);
  myTree->SetBranchAddress("ech_x", ech_x, &b_ech_x);
  myTree->SetBranchAddress("ech_y", ech_y, &b_ech_y);
  myTree->SetBranchAddress("ech_z", ech_z, &b_ech_z);
  myTree->SetBranchAddress("ec_m2", ec_m2, &b_ec_m2);
  myTree->SetBranchAddress("ec_m3", ec_m3, &b_ec_m3);
  myTree->SetBranchAddress("ec_m4", ec_m4, &b_ec_m4);
  myTree->SetBranchAddress("ec_c2", ec_c2, &b_ec_c2);
  myTree->SetBranchAddress("sc_part", &sc_part, &b_sc_part);
  myTree->SetBranchAddress("sc_sect", sc_sect, &b_sc_sect);
  myTree->SetBranchAddress("sc_hit", sc_hit, &b_sc_hit);
  myTree->SetBranchAddress("sc_pd", sc_pd, &b_sc_pd);
  myTree->SetBranchAddress("sc_stat", sc_stat, &b_sc_stat);
  myTree->SetBranchAddress("edep", edep, &b_edep);
  myTree->SetBranchAddress("sc_t", sc_t, &b_sc_t);
  myTree->SetBranchAddress("sc_r", sc_r, &b_sc_r);
  myTree->SetBranchAddress("sc_c2", sc_c2, &b_sc_c2);
  myTree->SetBranchAddress("cc_part", &cc_part, &b_cc_part);
  myTree->SetBranchAddress("cc_sect", cc_sect, &b_cc_sect);
  myTree->SetBranchAddress("cc_hit", cc_hit, &b_cc_hit);
  myTree->SetBranchAddress("cc_segm", cc_segm, &b_cc_segm);
  myTree->SetBranchAddress("nphe", nphe, &b_nphe);
  myTree->SetBranchAddress("cc_t", cc_t, &b_cc_t);
  myTree->SetBranchAddress("cc_r", cc_r, &b_cc_r);
  myTree->SetBranchAddress("cc_c2", cc_c2, &b_cc_c2);
  myTree->SetBranchAddress("lac_part", &lac_part, &b_lac_part);
  myTree->SetBranchAddress("lec_sect", &lec_sect, &b_lec_sect);
  myTree->SetBranchAddress("lec_hit", &lec_hit, &b_lec_hit);
  myTree->SetBranchAddress("lec_stat", &lec_stat, &b_lec_stat);
  myTree->SetBranchAddress("lec_etot", &lec_etot, &b_lec_etot);
  myTree->SetBranchAddress("lec_ein", &lec_ein, &b_lec_ein);
  myTree->SetBranchAddress("lec_t", &lec_t, &b_lec_t);
  myTree->SetBranchAddress("lec_r", &lec_r, &b_lec_r);
  myTree->SetBranchAddress("lec_x", &lec_x, &b_lec_x);
  myTree->SetBranchAddress("lec_y", &lec_y, &b_lec_y);
  myTree->SetBranchAddress("lec_z", &lec_z, &b_lec_z);
  myTree->SetBranchAddress("lec_c2", &lec_c2, &b_lec_c2);
  myTree->SetBranchAddress("st_part", &st_part, &b_st_part);
  myTree->SetBranchAddress("st_status", &st_status, &b_st_status);
  myTree->SetBranchAddress("st_time", &st_time, &b_st_time);
  myTree->SetBranchAddress("st_rtrk", &st_rtrk, &b_st_rtrk);

  myTree->SetBranchStatus("*", 1);
}
