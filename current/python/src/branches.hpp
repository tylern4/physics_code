#include "TTree.h"
#include "TChain.h"

/* My Branches */
Float_t W;
Float_t Q2;
Float_t MM;

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

void getBranches(TTree *myTree) {

  myTree->SetBranchAddress("npart", &npart);   // number of final particles
  myTree->SetBranchAddress("evntid", &evntid); // event number
  myTree->SetBranchAddress("evntclas", &evntclas);
  myTree->SetBranchAddress("q_l", &q_l);
  myTree->SetBranchAddress("t_l", &t_l);
  myTree->SetBranchAddress("tr_time", &tr_time);
  // number of particles in a single event (geometric particles)
  myTree->SetBranchAddress("gpart", &gpart);
  myTree->SetBranchAddress("id", &id); // particle ID of i'th element id[i]
  myTree->SetBranchAddress("stat", &stat);
  myTree->SetBranchAddress("dc", &dc);
  myTree->SetBranchAddress("cc", &cc);
  myTree->SetBranchAddress("sc", &sc);
  myTree->SetBranchAddress("ec", &ec);
  myTree->SetBranchAddress("lec", &lec);
  myTree->SetBranchAddress("p", &p); // momentum of i'th particle p[i] (GeV/C)
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
  myTree->SetBranchAddress("etot", &etot);
  myTree->SetBranchAddress("ec_ei", &ec_ei);
  myTree->SetBranchAddress("ec_eo", &ec_eo);
  myTree->SetBranchAddress("ec_t", &ec_t);
  myTree->SetBranchAddress("ec_r", &ec_r);
  myTree->SetBranchAddress("ech_x", &ech_x);
  myTree->SetBranchAddress("ech_y", &ech_y);
  myTree->SetBranchAddress("ech_z", &ech_z);
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

  myTree->SetBranchStatus("*", 1);
}
