/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#ifndef BRANCHES_H
#define BRANCHES_H
#include <iostream>
#include <vector>
#include "TChain.h"
#include "TTree.h"
#include "TTreeCache.h"
#include "TTreePerfStats.h"
#include "constants.hpp"

/*
******************************************************************************
*Tree    :h10       : h10
******************************************************************************
*Br    0 :npart     : npart/I
*Br    1 :evstat    : evstat/I
*Br    2 :intt      : intt/I
*Br    3 :evntid    : evntid/I
*Br    4 :evtype    : evtype/I
*Br    5 :evntclas  : evntclas/I
*Br    6 :evthel    : evthel/I
*Br    7 :evntclas2 : evntclas2/I
*Br    8 :q_l       : q_l/F
*Br    9 :t_l       : t_l/F
*Br   10 :tr_time   : tr_time/F
*Br   11 :rf_time1  : rf_time1/F
*Br   12 :rf_time2  : rf_time2/F
*Br   13 :gpart     : gpart/I
*Br   14 :id        : id[gpart]/I
*Br   15 :stat      : stat[gpart]/I
*Br   16 :dc        : dc[gpart]/I
*Br   17 :cc        : cc[gpart]/I
*Br   18 :sc        : sc[gpart]/I
*Br   19 :ec        : ec[gpart]/I
*Br   20 :lec       : lec[gpart]/I
*Br   21 :ccst      : ccst[gpart]/I
*Br   22 :p         : p[gpart]/F
*Br   23 :m         : m[gpart]/F
*Br   24 :q         : q[gpart]/I
*Br   25 :b         : b[gpart]/F
*Br   26 :cx        : cx[gpart]/F
*Br   27 :cy        : cy[gpart]/F
*Br   28 :cz        : cz[gpart]/F
*Br   29 :vx        : vx[gpart]/F
*Br   30 :vy        : vy[gpart]/F
*Br   31 :vz        : vz[gpart]/F
*Br   32 :dc_part   : dc_part/I
*Br   33 :dc_sect   : dc_sect[dc_part]/I
*Br   34 :dc_trk    : dc_trk[dc_part]/I
*Br   35 :dc_stat   : dc_stat[dc_part]/I
*Br   36 :dc_vx     : dc_vx[dc_part]/F
*Br   37 :dc_vy     : dc_vy[dc_part]/F
*Br   38 :dc_vz     : dc_vz[dc_part]/F
*Br   39 :dc_vr     : dc_vr[dc_part]/F
*Br   40 :dc_xsc    : dc_xsc[dc_part]/F
*Br   41 :dc_ysc    : dc_ysc[dc_part]/F
*Br   42 :dc_zsc    : dc_zsc[dc_part]/F
*Br   43 :dc_cxsc   : dc_cxsc[dc_part]/F
*Br   44 :dc_cysc   : dc_cysc[dc_part]/F
*Br   45 :dc_czsc   : dc_czsc[dc_part]/F
*Br   46 :dc_c2     : dc_c2[dc_part]/F
*Br   47 :ec_part   : ec_part/I
*Br   48 :ec_stat   : ec_stat[ec_part]/I
*Br   49 :ec_sect   : ec_sect[ec_part]/I
*Br   50 :ec_whol   : ec_whol[ec_part]/I
*Br   51 :ec_inst   : ec_inst[ec_part]/I
*Br   52 :ec_oust   : ec_oust[ec_part]/I
*Br   53 :etot      : etot[ec_part]/F
*Br   54 :ec_ei     : ec_ei[ec_part]/F
*Br   55 :ec_eo     : ec_eo[ec_part]/F
*Br   56 :ec_t      : ec_t[ec_part]/F
*Br   57 :ec_r      : ec_r[ec_part]/F
*Br   58 :ech_x     : ech_x[ec_part]/F
*Br   59 :ech_y     : ech_y[ec_part]/F
*Br   60 :ech_z     : ech_z[ec_part]/F
*Br   61 :ec_m2     : ec_m2[ec_part]/F
*Br   62 :ec_m3     : ec_m3[ec_part]/F
*Br   63 :ec_m4     : ec_m4[ec_part]/F
*Br   64 :ec_c2     : ec_c2[ec_part]/F
*Br   65 :sc_part   : sc_part/I
*Br   66 :sc_sect   : sc_sect[sc_part]/I
*Br   67 :sc_hit    : sc_hit[sc_part]/I
*Br   68 :sc_pd     : sc_pd[sc_part]/I
*Br   69 :sc_stat   : sc_stat[sc_part]/I
*Br   70 :edep      : edep[sc_part]/F
*Br   71 :sc_t      : sc_t[sc_part]/F
*Br   72 :sc_r      : sc_r[sc_part]/F
*Br   73 :sc_c2     : sc_c2[sc_part]/F
*Br   74 :cc_part   : cc_part/I
*Br   75 :cc_sect   : cc_sect[cc_part]/I
*Br   76 :cc_hit    : cc_hit[cc_part]/I
*Br   77 :cc_segm   : cc_segm[cc_part]/I
*Br   78 :nphe      : nphe[cc_part]/I
*Br   79 :cc_t      : cc_t[cc_part]/F
*Br   80 :cc_r      : cc_r[cc_part]/F
*Br   81 :cc_c2     : cc_c2[cc_part]/F
*Br   82 :lac_part  : lac_part/I
*Br   83 :lec_sect  : lec_sect[lac_part]/I
*Br   84 :lec_hit   : lec_hit[lac_part]/I
*Br   85 :lec_stat  : lec_stat[lac_part]/I
*Br   86 :lec_etot  : lec_etot[lac_part]/F
*Br   87 :lec_ein   : lec_ein[lac_part]/F
*Br   88 :lec_t     : lec_t[lac_part]/F
*Br   89 :lec_r     : lec_r[lac_part]/F
*Br   90 :lec_x     : lec_x[lac_part]/F
*Br   91 :lec_y     : lec_y[lac_part]/F
*Br   92 :lec_z     : lec_z[lac_part]/F
*Br   93 :lec_c2    : lec_c2[lac_part]/F
*Br   94 :st_part   : st_part/I
*Br   95 :st_status : st_status[st_part]/I
*Br   96 :st_time   : st_time[st_part]/F
*Br   97 :st_rtrk   : st_rtrk[st_part]/F
*/

class Branches {
 protected:
  std::shared_ptr<TTree> _tree;
  bool _MC = false;
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
  ////////////
  int _nprt;
  int _pidpart[MAX_PARTS];    //[nprt]
  float _xpart[MAX_PARTS];    //[nprt]
  float _ypart[MAX_PARTS];    //[nprt]
  float _zpart[MAX_PARTS];    //[nprt]
  float _epart[MAX_PARTS];    //[nprt]
  float _pxpart[MAX_PARTS];   //[nprt]
  float _pypart[MAX_PARTS];   //[nprt]
  float _pzpart[MAX_PARTS];   //[nprt]
  float _qpart[MAX_PARTS];    //[nprt]
  int _flagspart[MAX_PARTS];  //[nprt]

 public:
  Branches(std::shared_ptr<TChain> tree);
  Branches(std::shared_ptr<TChain> tree, bool MC);
  Branches(const Branches& b);
  ~Branches(){};
  void init();
  int npart();
  int evstat();
  int intt();
  int evntid();
  int evtype();
  int evntclas();
  int evthel();
  int evntclas2();
  int helicity();
  float q_l();
  float t_l();
  float tr_time();
  float rf_time1();
  float rf_time2();
  int gpart();
  int dc_part();
  int ec_part();
  int sc_part();
  int cc_part();
  int lac_part();
  int st_part();
  int nprt();

  std::vector<int> id();
  std::vector<int> stat();
  std::vector<int> dc();
  std::vector<int> cc();
  std::vector<int> sc();
  std::vector<int> ec();
  std::vector<int> lec();
  std::vector<int> ccst();
  std::vector<float> p();
  std::vector<float> px();
  std::vector<float> py();
  std::vector<float> pz();
  std::vector<float> m();
  std::vector<int> q();
  std::vector<float> b();
  std::vector<float> cx();
  std::vector<float> cy();
  std::vector<float> cz();
  std::vector<float> vx();
  std::vector<float> vy();
  std::vector<float> vz();
  std::vector<int> dc_sect();
  std::vector<int> dc_trk();
  std::vector<int> dc_stat();
  std::vector<float> dc_vx();
  std::vector<float> dc_vy();
  std::vector<float> dc_vz();
  std::vector<float> dc_vr();
  std::vector<float> dc_xsc();
  std::vector<float> dc_ysc();
  std::vector<float> dc_zsc();
  std::vector<float> dc_cxsc();
  std::vector<float> dc_cysc();
  std::vector<float> dc_czsc();
  std::vector<float> dc_c2();
  std::vector<int> ec_stat();
  std::vector<int> ec_sect();
  std::vector<int> ec_whol();
  std::vector<int> ec_inst();
  std::vector<int> ec_oust();
  std::vector<float> etot();
  std::vector<float> ec_ei();
  std::vector<float> ec_eo();
  std::vector<float> ec_t();
  std::vector<float> ec_r();
  std::vector<float> ech_x();
  std::vector<float> ech_y();
  std::vector<float> ech_z();
  std::vector<float> ec_m2();
  std::vector<float> ec_m3();
  std::vector<float> ec_m4();
  std::vector<float> ec_c2();
  std::vector<int> sc_sect();
  std::vector<int> sc_hit();
  std::vector<int> sc_pd();
  std::vector<int> sc_stat();
  std::vector<float> edep();
  std::vector<float> sc_t();
  std::vector<float> sc_r();
  std::vector<float> sc_c2();
  std::vector<int> cc_sect();
  std::vector<int> cc_hit();
  std::vector<int> cc_segm();
  std::vector<int> nphe();
  std::vector<float> cc_t();
  std::vector<float> cc_r();
  std::vector<float> cc_c2();
  //////////
  std::vector<int> pidpart();
  std::vector<float> xpart();
  std::vector<float> ypart();
  std::vector<float> zpart();
  std::vector<float> epart();
  std::vector<float> pxpart();
  std::vector<float> pypart();
  std::vector<float> pzpart();
  std::vector<float> qpart();

  int id(int i);
  int stat(int i);
  int dc(int i);
  int cc(int i);
  int sc(int i);
  int ec(int i);
  int lec(int i);
  int ccst(int i);
  float p(int i);
  float px(int i);
  float py(int i);
  float pz(int i);
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
  //////////
  int pidpart(int i);
  float xpart(int i);
  float ypart(int i);
  float zpart(int i);
  float epart(int i);
  float pxpart(int i);
  float pypart(int i);
  float pzpart(int i);
  float qpart(int i);
};

class BranchesMC : public Branches {
 public:
  BranchesMC(std::shared_ptr<TChain> tree);
};

#endif
