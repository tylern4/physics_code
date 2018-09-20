/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#ifndef BRANCHES_H
#define BRANCHES_H
#include <vector>
#include "TChain.h"
#include "constants.hpp"

class Branches {
 private:
  TChain* myTree;
  int _npart;
  int _evstat;
  int _intt;
  int _evntid;
  int _evtype;
  int _evntclas;
  int _evthel;
  int _evntclas2;
  float _q_l;
  float _t_l;
  float _tr_time;
  float _rf_time1;
  float _rf_time2;
  int _gpart;
  int _id[MAX_PARTS];    //[gpart]
  int _stat[MAX_PARTS];  //[gpart]
  int _dc[MAX_PARTS];    //[gpart]
  int _cc[MAX_PARTS];    //[gpart]
  int _sc[MAX_PARTS];    //[gpart]
  int _ec[MAX_PARTS];    //[gpart]
  int _lec[MAX_PARTS];   //[gpart]
  int _ccst[MAX_PARTS];  //[gpart]
  float _p[MAX_PARTS];   //[gpart]
  float _m[MAX_PARTS];   //[gpart]
  int _q[MAX_PARTS];     //[gpart]
  float _b[MAX_PARTS];   //[gpart]
  float _cx[MAX_PARTS];  //[gpart]
  float _cy[MAX_PARTS];  //[gpart]
  float _cz[MAX_PARTS];  //[gpart]
  float _vx[MAX_PARTS];  //[gpart]
  float _vy[MAX_PARTS];  //[gpart]
  float _vz[MAX_PARTS];  //[gpart]
  int _dc_part;
  int _dc_sect[MAX_PARTS];    //[dc_part]
  int _dc_trk[MAX_PARTS];     //[dc_part]
  int _dc_stat[MAX_PARTS];    //[dc_part]
  float _dc_vx[MAX_PARTS];    //[dc_part]
  float _dc_vy[MAX_PARTS];    //[dc_part]
  float _dc_vz[MAX_PARTS];    //[dc_part]
  float _dc_vr[MAX_PARTS];    //[dc_part]
  float _dc_xsc[MAX_PARTS];   //[dc_part]
  float _dc_ysc[MAX_PARTS];   //[dc_part]
  float _dc_zsc[MAX_PARTS];   //[dc_part]
  float _dc_cxsc[MAX_PARTS];  //[dc_part]
  float _dc_cysc[MAX_PARTS];  //[dc_part]
  float _dc_czsc[MAX_PARTS];  //[dc_part]
  float _dc_c2[MAX_PARTS];    //[dc_part]
  int _ec_part;
  int _ec_stat[MAX_PARTS];  //[ec_part]
  int _ec_sect[MAX_PARTS];  //[ec_part]
  int _ec_whol[MAX_PARTS];  //[ec_part]
  int _ec_inst[MAX_PARTS];  //[ec_part]
  int _ec_oust[MAX_PARTS];  //[ec_part]
  float _etot[MAX_PARTS];   //[ec_part]
  float _ec_ei[MAX_PARTS];  //[ec_part]
  float _ec_eo[MAX_PARTS];  //[ec_part]
  float _ec_t[MAX_PARTS];   //[ec_part]
  float _ec_r[MAX_PARTS];   //[ec_part]
  float _ech_x[MAX_PARTS];  //[ec_part]
  float _ech_y[MAX_PARTS];  //[ec_part]
  float _ech_z[MAX_PARTS];  //[ec_part]
  float _ec_m2[MAX_PARTS];  //[ec_part]
  float _ec_m3[MAX_PARTS];  //[ec_part]
  float _ec_m4[MAX_PARTS];  //[ec_part]
  float _ec_c2[MAX_PARTS];  //[ec_part]
  int _sc_part;
  int _sc_sect[MAX_PARTS];  //[sc_part]
  int _sc_hit[MAX_PARTS];   //[sc_part]
  int _sc_pd[MAX_PARTS];    //[sc_part]
  int _sc_stat[MAX_PARTS];  //[sc_part]
  float _edep[MAX_PARTS];   //[sc_part]
  float _sc_t[MAX_PARTS];   //[sc_part]
  float _sc_r[MAX_PARTS];   //[sc_part]
  float _sc_c2[MAX_PARTS];  //[sc_part]
  int _cc_part;
  int _cc_sect[MAX_PARTS];  //[cc_part]
  int _cc_hit[MAX_PARTS];   //[cc_part]
  int _cc_segm[MAX_PARTS];  //[cc_part]
  int _nphe[MAX_PARTS];     //[cc_part]
  float _cc_t[MAX_PARTS];   //[cc_part]
  float _cc_r[MAX_PARTS];   //[cc_part]
  float _cc_c2[MAX_PARTS];  //[cc_part]
  int _lac_part;
  int _lec_sect[MAX_PARTS];    //[lac_part]
  int _lec_hit[MAX_PARTS];     //[lac_part]
  int _lec_stat[MAX_PARTS];    //[lac_part]
  float _lec_etot[MAX_PARTS];  //[lac_part]
  float _lec_ein[MAX_PARTS];   //[lac_part]
  float _lec_t[MAX_PARTS];     //[lac_part]
  float _lec_r[MAX_PARTS];     //[lac_part]
  float _lec_x[MAX_PARTS];     //[lac_part]
  float _lec_y[MAX_PARTS];     //[lac_part]
  float _lec_z[MAX_PARTS];     //[lac_part]
  float _lec_c2[MAX_PARTS];    //[lac_part]
  int _st_part;
  int _st_status[MAX_PARTS];  //[st_part]
  float _st_time[MAX_PARTS];  //[st_part]
  float _st_rtrk[MAX_PARTS];  //[st_part]
 public:
  Branches(TChain* tree);
  Branches(const Branches& b);
  ~Branches();
  void init();
  int npart();
  int evstat();
  int intt();
  int evntid();
  int evtype();
  int evntclas();
  int evthel();
  int evntclas2();
  float q_l();
  float t_l();
  float tr_time();
  float rf_time1();
  float rf_time2();
  int gpart();
  std::vector<int> id();
  std::vector<int> stat();
  std::vector<int> dc();
  std::vector<int> cc();
  std::vector<int> sc();
  std::vector<int> ec();
  std::vector<int> lec();
  std::vector<int> ccst();
  std::vector<float> p();
  std::vector<float> m();
  std::vector<int> q();
  std::vector<float> b();
  std::vector<float> cx();
  std::vector<float> cy();
  std::vector<float> cz();
  std::vector<float> vx();
  std::vector<float> vy();
  std::vector<float> vz();

  int dc_part();
  int ec_part();
  int sc_part();
  int cc_part();
  int lac_part();
  int st_part();

  int id(int i);
  int stat(int i);
  int dc(int i);
  int cc(int i);
  int sc(int i);
  int ec(int i);
  int lec(int i);
  int ccst(int i);
  float p(int i);
  float m(int i);
  int q(int i);
  float b(int i);
  float cx(int i);
  float cy(int i);
  float cz(int i);
  float vx(int i);
  float vy(int i);
  float vz(int i);
  int dc_sect(int i);
  int dc_trk(int i);
  int dc_stat(int i);
  float dc_vx(int i);
  float dc_vy(int i);
  float dc_vz(int i);
  float dc_vr(int i);
  float dc_xsc(int i);
  float dc_ysc(int i);
  float dc_zsc(int i);
  float dc_cxsc(int i);
  float dc_cysc(int i);
  float dc_czsc(int i);
  float dc_c2(int i);
  int ec_stat(int i);
  int ec_sect(int i);
  int ec_whol(int i);
  int ec_inst(int i);
  int ec_oust(int i);
  float etot(int i);
  float ec_ei(int i);
  float ec_eo(int i);
  float ec_t(int i);
  float ec_r(int i);
  float ech_x(int i);
  float ech_y(int i);
  float ech_z(int i);
  float ec_m2(int i);
  float ec_m3(int i);
  float ec_m4(int i);
  float ec_c2(int i);

  int sc_sect(int i);
  int sc_hit(int i);
  int sc_pd(int i);
  int sc_stat(int i);
  float edep(int i);
  float sc_t(int i);
  float sc_r(int i);
  float sc_c2(int i);

  int cc_sect(int i);
  int cc_hit(int i);
  int cc_segm(int i);
  int nphe(int i);
  float cc_t(int i);
  float cc_r(int i);
  float cc_c2(int i);
};

#endif
