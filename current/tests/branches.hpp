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

extern UChar_t npart;
// extern UChar_t evstat;
extern UInt_t evntid;
extern Char_t evtype;
extern Char_t evntclas;
extern Char_t evthel;
extern Int_t evntclas2;
extern Float_t q_l;
extern Float_t t_l;
extern Float_t tr_time;
extern Float_t rf_time1;
extern Float_t rf_time2;

extern Int_t gpart;
extern Short_t id[MAX_PARTS];  //[gpart]
extern Char_t stat[MAX_PARTS]; //[gpart]
extern UChar_t dc[MAX_PARTS];  //[gpart]
extern UChar_t cc[MAX_PARTS];  //[gpart]
extern UChar_t sc[MAX_PARTS];  //[gpart]
extern UChar_t ec[MAX_PARTS];  //[gpart]
extern UChar_t lec[MAX_PARTS]; //[gpart]
extern Float_t p[MAX_PARTS];   //[gpart]
extern Float_t m[MAX_PARTS];   //[gpart]
extern Char_t q[MAX_PARTS];    //[gpart]
extern Float_t b[MAX_PARTS];   //[gpart]
extern Float_t cx[MAX_PARTS];  //[gpart]
extern Float_t cy[MAX_PARTS];  //[gpart]
extern Float_t cz[MAX_PARTS];  //[gpart]
extern Float_t vx[MAX_PARTS];  //[gpart]
extern Float_t vy[MAX_PARTS];  //[gpart]
extern Float_t vz[MAX_PARTS];  //[gpart]

// Each of the folowing has multiple parts

extern Int_t dc_part;
extern UChar_t dc_sect[MAX_PARTS]; //[dc_part]
extern UChar_t dc_trk[MAX_PARTS];  //[dc_part]
extern Char_t dc_stat[MAX_PARTS];  //[dc_part]
extern Float_t dc_xsc[MAX_PARTS];  //[dc_part]
extern Float_t dc_ysc[MAX_PARTS];  //[dc_part]
extern Float_t dc_zsc[MAX_PARTS];  //[dc_part]
extern Float_t dc_cxsc[MAX_PARTS]; //[dc_part]
extern Float_t dc_cysc[MAX_PARTS]; //[dc_part]
extern Float_t dc_czsc[MAX_PARTS]; //[dc_part]
extern Float_t dc_xec[MAX_PARTS];  //[dc_part]
extern Float_t dc_yec[MAX_PARTS];  //[dc_part]
extern Float_t dc_zec[MAX_PARTS];  //[dc_part]
extern Float_t dc_thcc[MAX_PARTS]; //[dc_part]
extern Float_t dc_c2[MAX_PARTS];   //[dc_part]

extern Int_t ec_part;
extern UShort_t ec_stat[MAX_PARTS]; //[ec_part]
extern UChar_t ec_sect[MAX_PARTS];  //[ec_part]
extern Int_t ec_whol[MAX_PARTS];    //[ec_part]
// extern Int_t   ec_inst[MAX_PARTS];   //[ec_part]
// extern Int_t   ec_oust[MAX_PARTS];   //[ec_part]
extern Float_t etot[MAX_PARTS];  //[ec_part]
extern Float_t ec_ei[MAX_PARTS]; //[ec_part]
extern Float_t ec_eo[MAX_PARTS]; //[ec_part]
extern Float_t ec_t[MAX_PARTS];  //[ec_part]
extern Float_t ec_r[MAX_PARTS];  //[ec_part]
extern Float_t ech_x[MAX_PARTS]; //[ec_part]
extern Float_t ech_y[MAX_PARTS]; //[ec_part]
extern Float_t ech_z[MAX_PARTS]; //[ec_part]
// extern Float_t ec_m2[MAX_PARTS];   //[ec_part]
// extern Float_t ec_m3[MAX_PARTS];   //[ec_part]
// extern Float_t ec_m4[MAX_PARTS];   //[ec_part]
extern Float_t ec_c2[MAX_PARTS]; //[ec_part]

extern Int_t sc_part;
extern UChar_t sc_sect[MAX_PARTS]; //[sc_part]
extern UChar_t sc_hit[MAX_PARTS];  //[sc_part]
extern UChar_t sc_pd[MAX_PARTS];   //[sc_part]
extern UChar_t sc_stat[MAX_PARTS]; //[sc_part]
extern Float_t edep[MAX_PARTS];    //[sc_part]
extern Float_t sc_t[MAX_PARTS];    //[sc_part]
extern Float_t sc_r[MAX_PARTS];    //[sc_part]
extern Float_t sc_c2[MAX_PARTS];   //[sc_part]

extern Int_t cc_part;
extern UChar_t cc_sect[MAX_PARTS]; //[cc_part]
extern UChar_t cc_hit[MAX_PARTS];  //[cc_part]
extern Int_t cc_segm[MAX_PARTS];   //[cc_part]
extern UShort_t nphe[MAX_PARTS];   //[cc_part]
extern Float_t cc_t[MAX_PARTS];    //[cc_part]
extern Float_t cc_r[MAX_PARTS];    //[cc_part]
extern Float_t cc_c2[MAX_PARTS];   //[cc_part]

/*extern Int_t   lac_part;
extern Int_t   lec_sect[MAX_PARTS];   //[lac_part]
extern Int_t   lec_hit[MAX_PARTS];   //[lac_part]
extern Int_t   lec_stat[MAX_PARTS];   //[lac_part]
extern Float_t lec_etot[MAX_PARTS];   //[lac_part]
extern Float_t lec_t[MAX_PARTS];   //[lac_part]
extern Float_t lec_r[MAX_PARTS];   //[lac_part]
extern Float_t lec_x[MAX_PARTS];   //[lac_part]
extern Float_t lec_y[MAX_PARTS];   //[lac_part]
extern Float_t lec_z[MAX_PARTS];   //[lac_part]
extern Float_t lec_c2[MAX_PARTS];   //[lac_part]
*/

void getMorebranchs(TTree *myTree);
void getBranches(TTree *myTree);

#endif
