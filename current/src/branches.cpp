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

UChar_t npart;
// UChar_t evstat;
UInt_t evntid;
Char_t evtype;
Char_t evntclas;
Char_t evthel;
Int_t evntclas2;
Float_t q_l;
Float_t t_l;
Float_t tr_time;
Float_t rf_time1;
Float_t rf_time2;

Int_t gpart;
Short_t id[MAX_PARTS];  //[gpart]
Char_t stat[MAX_PARTS]; //[gpart]
UChar_t dc[MAX_PARTS];  //[gpart]
UChar_t cc[MAX_PARTS];  //[gpart]
UChar_t sc[MAX_PARTS];  //[gpart]
UChar_t ec[MAX_PARTS];  //[gpart]
UChar_t lec[MAX_PARTS]; //[gpart]
Float_t p[MAX_PARTS];   //[gpart]
Float_t m[MAX_PARTS];   //[gpart]
Char_t q[MAX_PARTS];    //[gpart]
Float_t b[MAX_PARTS];   //[gpart]
Float_t cx[MAX_PARTS];  //[gpart]
Float_t cy[MAX_PARTS];  //[gpart]
Float_t cz[MAX_PARTS];  //[gpart]
Float_t vx[MAX_PARTS];  //[gpart]
Float_t vy[MAX_PARTS];  //[gpart]
Float_t vz[MAX_PARTS];  //[gpart]

// Each of the folowing has multiple parts

Int_t dc_part;
UChar_t dc_sect[MAX_PARTS]; //[dc_part]
UChar_t dc_trk[MAX_PARTS];  //[dc_part]
Char_t dc_stat[MAX_PARTS];  //[dc_part]
Float_t dc_xsc[MAX_PARTS];  //[dc_part]
Float_t dc_ysc[MAX_PARTS];  //[dc_part]
Float_t dc_zsc[MAX_PARTS];  //[dc_part]
Float_t dc_cxsc[MAX_PARTS]; //[dc_part]
Float_t dc_cysc[MAX_PARTS]; //[dc_part]
Float_t dc_czsc[MAX_PARTS]; //[dc_part]
Float_t dc_xec[MAX_PARTS];  //[dc_part]
Float_t dc_yec[MAX_PARTS];  //[dc_part]
Float_t dc_zec[MAX_PARTS];  //[dc_part]
Float_t dc_thcc[MAX_PARTS]; //[dc_part]
Float_t dc_c2[MAX_PARTS];   //[dc_part]

Int_t ec_part;
UShort_t ec_stat[MAX_PARTS]; //[ec_part]
UChar_t ec_sect[MAX_PARTS];  //[ec_part]
Int_t ec_whol[MAX_PARTS];    //[ec_part]
// Int_t   ec_inst[MAX_PARTS];   //[ec_part]
// Int_t   ec_oust[MAX_PARTS];   //[ec_part]
Float_t etot[MAX_PARTS];  //[ec_part]
Float_t ec_ei[MAX_PARTS]; //[ec_part]
Float_t ec_eo[MAX_PARTS]; //[ec_part]
Float_t ec_t[MAX_PARTS];  //[ec_part]
Float_t ec_r[MAX_PARTS];  //[ec_part]
Float_t ech_x[MAX_PARTS]; //[ec_part]
Float_t ech_y[MAX_PARTS]; //[ec_part]
Float_t ech_z[MAX_PARTS]; //[ec_part]
// Float_t ec_m2[MAX_PARTS];   //[ec_part]
// Float_t ec_m3[MAX_PARTS];   //[ec_part]
// Float_t ec_m4[MAX_PARTS];   //[ec_part]
Float_t ec_c2[MAX_PARTS]; //[ec_part]

Int_t sc_part;
UChar_t sc_sect[MAX_PARTS]; //[sc_part]
UChar_t sc_hit[MAX_PARTS];  //[sc_part]
UChar_t sc_pd[MAX_PARTS];   //[sc_part]
UChar_t sc_stat[MAX_PARTS]; //[sc_part]
Float_t edep[MAX_PARTS];    //[sc_part]
Float_t sc_t[MAX_PARTS];    //[sc_part]
Float_t sc_r[MAX_PARTS];    //[sc_part]
Float_t sc_c2[MAX_PARTS];   //[sc_part]

Int_t cc_part;
UChar_t cc_sect[MAX_PARTS]; //[cc_part]
UChar_t cc_hit[MAX_PARTS];  //[cc_part]
Int_t cc_segm[MAX_PARTS];   //[cc_part]
UShort_t nphe[MAX_PARTS];   //[cc_part]
Float_t cc_t[MAX_PARTS];    //[cc_part]
Float_t cc_r[MAX_PARTS];    //[cc_part]
Float_t cc_c2[MAX_PARTS];   //[cc_part]

/*Int_t   lac_part;
Int_t   lec_sect[MAX_PARTS];   //[lac_part]
Int_t   lec_hit[MAX_PARTS];   //[lac_part]
Int_t   lec_stat[MAX_PARTS];   //[lac_part]
Float_t lec_etot[MAX_PARTS];   //[lac_part]
Float_t lec_t[MAX_PARTS];   //[lac_part]
Float_t lec_r[MAX_PARTS];   //[lac_part]
Float_t lec_x[MAX_PARTS];   //[lac_part]
Float_t lec_y[MAX_PARTS];   //[lac_part]
Float_t lec_z[MAX_PARTS];   //[lac_part]
Float_t lec_c2[MAX_PARTS];   //[lac_part]
*/

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

  myTree->SetBranchAddress("npart", &npart); // number of final particles
  // myTree->SetBranchAddress("evstat", &evstat);
  myTree->SetBranchAddress("evntid", &evntid); // event number
  myTree->SetBranchAddress("evntclas", &evntclas);
  myTree->SetBranchAddress("q_l", &q_l);
  myTree->SetBranchAddress("t_l", &t_l);
  myTree->SetBranchAddress("tr_time", &tr_time);
  myTree->SetBranchAddress(
      "gpart",
      &gpart); // number of particles in a single event (geometric particles)
  myTree->SetBranchAddress("id", &id); // particle ID of i'th element id[i]
  myTree->SetBranchAddress("stat", &stat);
  myTree->SetBranchAddress("dc", &dc);
  myTree->SetBranchAddress("cc", &cc);
  myTree->SetBranchAddress("sc", &sc);
  myTree->SetBranchAddress("ec", &ec);
  myTree->SetBranchAddress("lec", &lec);
  myTree->SetBranchAddress("p", &p); // momentum of i'th particle p[i] (GeV/C)
  // mass of the particle m[i] (GeV/C)
  // myTree->SetBranchAddress("m", &m);
  // charge of i'th particle q[i] (charge in e's 1,0,-1)
  myTree->SetBranchAddress("q", &q);
  // Velocity of i'th particle b[i] (in terms of c) ie. Beta
  myTree->SetBranchAddress("b", &b);
  myTree->SetBranchAddress("cx", &cx); // X direction cosine at origin
  myTree->SetBranchAddress("cy", &cy); // Y direction cosine at origin
  myTree->SetBranchAddress("cz", &cz); // Z direction cosine at origin
  myTree->SetBranchAddress("vx", &vx); // X coordinate of vertex (cm)
  myTree->SetBranchAddress("vy", &vy); // y coordinate of vertex (cm)
  myTree->SetBranchAddress("vz", &vz); // z coordinate of vertex (cm)
  myTree->SetBranchAddress("dc_part", &dc_part);
  myTree->SetBranchAddress("dc_sect", &dc_sect);
  myTree->SetBranchAddress("dc_trk", &dc_trk);
  myTree->SetBranchAddress("dc_stat", &dc_stat);
  myTree->SetBranchAddress("dc_xsc", &dc_xsc);
  myTree->SetBranchAddress("dc_ysc", &dc_ysc);
  myTree->SetBranchAddress("dc_zsc", &dc_zsc);
  myTree->SetBranchAddress("dc_cxsc", &dc_cxsc);
  myTree->SetBranchAddress("dc_cysc", &dc_cysc);
  myTree->SetBranchAddress("dc_czsc", &dc_czsc);
  myTree->SetBranchAddress("dc_c2", &dc_c2);
  myTree->SetBranchAddress("ec_part", &ec_part);
  myTree->SetBranchAddress("ec_stat", &ec_stat);
  myTree->SetBranchAddress("ec_sect", &ec_sect);
  myTree->SetBranchAddress("ec_whol", &ec_whol);
  // myTree->SetBranchAddress("ec_inst", &ec_inst);
  // myTree->SetBranchAddress("ec_oust", &ec_oust);
  myTree->SetBranchAddress("etot", &etot);
  myTree->SetBranchAddress("ec_ei", &ec_ei);
  myTree->SetBranchAddress("ec_eo", &ec_eo);
  myTree->SetBranchAddress("ec_t", &ec_t);
  myTree->SetBranchAddress("ec_r", &ec_r);
  myTree->SetBranchAddress("ech_x", &ech_x);
  myTree->SetBranchAddress("ech_y", &ech_y);
  myTree->SetBranchAddress("ech_z", &ech_z);
  // myTree->SetBranchAddress("ec_m2", &ec_m2);
  // myTree->SetBranchAddress("ec_m3", &ec_m3);
  // myTree->SetBranchAddress("ec_m4", &ec_m4);
  myTree->SetBranchAddress("ec_c2", &ec_c2);
  myTree->SetBranchAddress("sc_part", &sc_part);
  myTree->SetBranchAddress("sc_sect", &sc_sect);
  myTree->SetBranchAddress("sc_hit", &sc_hit);
  myTree->SetBranchAddress("sc_pd", &sc_pd);
  myTree->SetBranchAddress("sc_stat", &sc_stat);
  myTree->SetBranchAddress("edep", &edep);
  myTree->SetBranchAddress("sc_t", &sc_t);
  myTree->SetBranchAddress("sc_r", &sc_r);
  myTree->SetBranchAddress("sc_c2", &sc_c2);
  myTree->SetBranchAddress("cc_part", &cc_part);
  myTree->SetBranchAddress("cc_sect", &cc_sect);
  myTree->SetBranchAddress("cc_hit", &cc_hit);
  myTree->SetBranchAddress("cc_segm", &cc_segm);
  myTree->SetBranchAddress("nphe", &nphe);
  myTree->SetBranchAddress("cc_t", &cc_t);
  myTree->SetBranchAddress("cc_r", &cc_r);
  myTree->SetBranchAddress("cc_c2", &cc_c2);
  /*myTree->SetBranchAddress("lac_part", &lac_part);
  myTree->SetBranchAddress("lec_sect", &lec_sect);
  myTree->SetBranchAddress("lec_hit", &lec_hit);
  myTree->SetBranchAddress("lec_stat", &lec_stat);
  myTree->SetBranchAddress("lec_etot", &lec_etot);
  myTree->SetBranchAddress("lec_t", &lec_t);
  myTree->SetBranchAddress("lec_r", &lec_r);
  myTree->SetBranchAddress("lec_x", &lec_x);
  myTree->SetBranchAddress("lec_y", &lec_y);
  myTree->SetBranchAddress("lec_z", &lec_z);
  myTree->SetBranchAddress("lec_c2", &lec_c2); */

  myTree->SetBranchStatus("*", 1);
}
