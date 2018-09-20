/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#include "branches.hpp"

Branches::Branches(TChain* tree) {
  myTree = tree;
  init();
}

Branches::Branches(const Branches& b) {
  myTree = b.myTree;
  init();
}

void Branches::init() {
  myTree->SetBranchAddress("npart", &_npart);
  myTree->SetBranchAddress("evstat", &_evstat);
  myTree->SetBranchAddress("intt", &_intt);
  myTree->SetBranchAddress("evntid", &_evntid);
  myTree->SetBranchAddress("evtype", &_evtype);
  myTree->SetBranchAddress("evntclas", &_evntclas);
  myTree->SetBranchAddress("evthel", &_evthel);
  myTree->SetBranchAddress("evntclas2", &_evntclas2);
  myTree->SetBranchAddress("q_l", &_q_l);
  myTree->SetBranchAddress("t_l", &_t_l);
  myTree->SetBranchAddress("tr_time", &_tr_time);
  myTree->SetBranchAddress("rf_time1", &_rf_time1);
  myTree->SetBranchAddress("rf_time2", &_rf_time2);
  myTree->SetBranchAddress("gpart", &_gpart);
  myTree->SetBranchAddress("id", _id);
  myTree->SetBranchAddress("stat", _stat);
  myTree->SetBranchAddress("dc", _dc);
  myTree->SetBranchAddress("cc", _cc);
  myTree->SetBranchAddress("sc", _sc);
  myTree->SetBranchAddress("ec", _ec);
  myTree->SetBranchAddress("lec", _lec);
  myTree->SetBranchAddress("ccst", _ccst);
  myTree->SetBranchAddress("p", _p);
  myTree->SetBranchAddress("m", _m);
  myTree->SetBranchAddress("q", _q);
  myTree->SetBranchAddress("b", _b);
  myTree->SetBranchAddress("cx", _cx);
  myTree->SetBranchAddress("cy", _cy);
  myTree->SetBranchAddress("cz", _cz);
  myTree->SetBranchAddress("vx", _vx);
  myTree->SetBranchAddress("vy", _vy);
  myTree->SetBranchAddress("vz", _vz);
  myTree->SetBranchAddress("dc_part", &_dc_part);
  myTree->SetBranchAddress("dc_sect", _dc_sect);
  myTree->SetBranchAddress("dc_trk", _dc_trk);
  myTree->SetBranchAddress("dc_stat", _dc_stat);
  myTree->SetBranchAddress("dc_vx", _dc_vx);
  myTree->SetBranchAddress("dc_vy", _dc_vy);
  myTree->SetBranchAddress("dc_vz", _dc_vz);
  myTree->SetBranchAddress("dc_vr", _dc_vr);
  myTree->SetBranchAddress("dc_xsc", _dc_xsc);
  myTree->SetBranchAddress("dc_ysc", _dc_ysc);
  myTree->SetBranchAddress("dc_zsc", _dc_zsc);
  myTree->SetBranchAddress("dc_cxsc", _dc_cxsc);
  myTree->SetBranchAddress("dc_cysc", _dc_cysc);
  myTree->SetBranchAddress("dc_czsc", _dc_czsc);
  myTree->SetBranchAddress("dc_c2", _dc_c2);
  myTree->SetBranchAddress("ec_part", &_ec_part);
  myTree->SetBranchAddress("ec_stat", _ec_stat);
  myTree->SetBranchAddress("ec_sect", _ec_sect);
  myTree->SetBranchAddress("ec_whol", _ec_whol);
  myTree->SetBranchAddress("ec_inst", _ec_inst);
  myTree->SetBranchAddress("ec_oust", _ec_oust);
  myTree->SetBranchAddress("etot", _etot);
  myTree->SetBranchAddress("ec_ei", _ec_ei);
  myTree->SetBranchAddress("ec_eo", _ec_eo);
  myTree->SetBranchAddress("ec_t", _ec_t);
  myTree->SetBranchAddress("ec_r", _ec_r);
  myTree->SetBranchAddress("ech_x", _ech_x);
  myTree->SetBranchAddress("ech_y", _ech_y);
  myTree->SetBranchAddress("ech_z", _ech_z);
  myTree->SetBranchAddress("ec_m2", _ec_m2);
  myTree->SetBranchAddress("ec_m3", _ec_m3);
  myTree->SetBranchAddress("ec_m4", _ec_m4);
  myTree->SetBranchAddress("ec_c2", _ec_c2);
  myTree->SetBranchAddress("sc_part", &_sc_part);
  myTree->SetBranchAddress("sc_sect", _sc_sect);
  myTree->SetBranchAddress("sc_hit", _sc_hit);
  myTree->SetBranchAddress("sc_pd", _sc_pd);
  myTree->SetBranchAddress("sc_stat", _sc_stat);
  myTree->SetBranchAddress("edep", _edep);
  myTree->SetBranchAddress("sc_t", _sc_t);
  myTree->SetBranchAddress("sc_r", _sc_r);
  myTree->SetBranchAddress("sc_c2", _sc_c2);
  myTree->SetBranchAddress("cc_part", &_cc_part);
  myTree->SetBranchAddress("cc_sect", _cc_sect);
  myTree->SetBranchAddress("cc_hit", _cc_hit);
  myTree->SetBranchAddress("cc_segm", _cc_segm);
  myTree->SetBranchAddress("nphe", _nphe);
  myTree->SetBranchAddress("cc_t", _cc_t);
  myTree->SetBranchAddress("cc_r", _cc_r);
  myTree->SetBranchAddress("cc_c2", _cc_c2);
  myTree->SetBranchAddress("lac_part", &_lac_part);
  myTree->SetBranchAddress("lec_sect", &_lec_sect);
  myTree->SetBranchAddress("lec_hit", &_lec_hit);
  myTree->SetBranchAddress("lec_stat", &_lec_stat);
  myTree->SetBranchAddress("lec_etot", &_lec_etot);
  myTree->SetBranchAddress("lec_ein", &_lec_ein);
  myTree->SetBranchAddress("lec_t", &_lec_t);
  myTree->SetBranchAddress("lec_r", &_lec_r);
  myTree->SetBranchAddress("lec_x", &_lec_x);
  myTree->SetBranchAddress("lec_y", &_lec_y);
  myTree->SetBranchAddress("lec_z", &_lec_z);
  myTree->SetBranchAddress("lec_c2", &_lec_c2);
  myTree->SetBranchAddress("st_part", &_st_part);
  myTree->SetBranchAddress("st_status", &_st_status);
  myTree->SetBranchAddress("st_time", &_st_time);
  myTree->SetBranchAddress("st_rtrk", &_st_rtrk);

  myTree->SetBranchStatus("*", 1);
}

int Branches::npart() { return _npart; }
int Branches::evstat() { return _evstat; }
int Branches::intt() { return _intt; }
int Branches::evntid() { return _evntid; }
int Branches::evtype() { return _evtype; }
int Branches::evntclas() { return _evntclas; }
int Branches::evthel() { return _evthel; }
int Branches::evntclas2() { return _evntclas2; }
float Branches::q_l() { return _q_l; }
float Branches::t_l() { return _t_l; }
float Branches::tr_time() { return _tr_time; }
float Branches::rf_time1() { return _rf_time1; }
float Branches::rf_time2() { return _rf_time2; }
int Branches::gpart() { return _gpart; }

int Branches::dc_part() { return _dc_part; }
int Branches::ec_part() { return _ec_part; }
int Branches::sc_part() { return _sc_part; }
int Branches::cc_part() { return _cc_part; }
int Branches::lac_part() { return _lac_part; }
int Branches::st_part() { return _st_part; }

int Branches::id(int i) {
  try {
    return _id[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  // [gpart]
int Branches::stat(int i) {
  try {
    return _stat[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  // [gpart]
int Branches::dc(int i) {
  try {
    return _dc[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  // [gpart]
int Branches::cc(int i) {
  try {
    return _cc[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  // [gpart]
int Branches::sc(int i) {
  try {
    return _sc[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  // [gpart]
int Branches::ec(int i) {
  try {
    return _ec[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  // [gpart]
int Branches::lec(int i) {
  try {
    return _lec[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  // [gpart]
int Branches::ccst(int i) {
  try {
    return _ccst[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  // [gpart]
float Branches::p(int i) {
  try {
    return _p[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  // [gpart]
float Branches::m(int i) {
  try {
    return _m[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  // [gpart]
int Branches::q(int i) {
  try {
    return _q[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  // [gpart]
float Branches::b(int i) {
  try {
    return _b[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  // [gpart]
float Branches::cx(int i) {
  try {
    return _cx[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  // [gpart]
float Branches::cy(int i) {
  try {
    return _cy[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  // [gpart]
float Branches::cz(int i) {
  try {
    return _cz[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  // [gpart]
float Branches::vx(int i) {
  try {
    return _vx[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  // [gpart]
float Branches::vy(int i) {
  try {
    return _vy[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  // [gpart]
float Branches::vz(int i) {
  try {
    return _vz[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  // [gpart]

int Branches::dc_sect(int i) {
  try {
    return _dc_sect[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[dc_part]
int Branches::dc_trk(int i) {
  try {
    return _dc_trk[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[dc_part]
int Branches::dc_stat(int i) {
  try {
    return _dc_stat[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[dc_part]
float Branches::dc_vx(int i) {
  try {
    return _dc_vx[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[dc_part]
float Branches::dc_vy(int i) {
  try {
    return _dc_vy[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[dc_part]
float Branches::dc_vz(int i) {
  try {
    return _dc_vz[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[dc_part]
float Branches::dc_vr(int i) {
  try {
    return _dc_vr[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[dc_part]
float Branches::dc_xsc(int i) {
  try {
    return _dc_xsc[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[dc_part]
float Branches::dc_ysc(int i) {
  try {
    return _dc_ysc[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[dc_part]
float Branches::dc_zsc(int i) {
  try {
    return _dc_zsc[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[dc_part]
float Branches::dc_cxsc(int i) {
  try {
    return _dc_cxsc[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[dc_part]
float Branches::dc_cysc(int i) {
  try {
    return _dc_cysc[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[dc_part]
float Branches::dc_czsc(int i) {
  try {
    return _dc_czsc[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[dc_part]
float Branches::dc_c2(int i) {
  try {
    return _dc_c2[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[dc_part]

int Branches::ec_stat(int i) {
  try {
    return _ec_stat[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[ec_part]
int Branches::ec_sect(int i) {
  try {
    return _ec_sect[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[ec_part]
int Branches::ec_whol(int i) {
  try {
    return _ec_whol[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[ec_part]
int Branches::ec_inst(int i) {
  try {
    return _ec_inst[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[ec_part]
int Branches::ec_oust(int i) {
  try {
    return _ec_oust[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[ec_part]
float Branches::etot(int i) {
  try {
    return _etot[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[ec_part]
float Branches::ec_ei(int i) {
  try {
    return _ec_ei[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[ec_part]
float Branches::ec_eo(int i) {
  try {
    return _ec_eo[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[ec_part]
float Branches::ec_t(int i) {
  try {
    return _ec_t[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[ec_part]
float Branches::ec_r(int i) {
  try {
    return _ec_r[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[ec_part]
float Branches::ech_x(int i) {
  try {
    return _ech_x[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[ec_part]
float Branches::ech_y(int i) {
  try {
    return _ech_y[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[ec_part]
float Branches::ech_z(int i) {
  try {
    return _ech_z[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[ec_part]
float Branches::ec_m2(int i) {
  try {
    return _ec_m2[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[ec_part]
float Branches::ec_m3(int i) {
  try {
    return _ec_m3[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[ec_part]
float Branches::ec_m4(int i) {
  try {
    return _ec_m4[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[ec_part]
float Branches::ec_c2(int i) {
  try {
    return _ec_c2[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[ec_part]

int Branches::sc_sect(int i) {
  try {
    return _sc_sect[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[sc_part]
int Branches::sc_hit(int i) {
  try {
    return _sc_hit[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[sc_part]
int Branches::sc_pd(int i) {
  try {
    return _sc_pd[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[sc_part]
int Branches::sc_stat(int i) {
  try {
    return _sc_stat[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[sc_part]
float Branches::edep(int i) {
  try {
    return _edep[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[sc_part]
float Branches::sc_t(int i) {
  try {
    return _sc_t[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[sc_part]
float Branches::sc_r(int i) {
  try {
    return _sc_r[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[sc_part]
float Branches::sc_c2(int i) {
  try {
    return _sc_c2[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[sc_part]

int Branches::cc_sect(int i) {
  try {
    return _cc_sect[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[cc_part]
int Branches::cc_hit(int i) {
  try {
    return _cc_hit[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[cc_part]
int Branches::cc_segm(int i) {
  try {
    return _cc_segm[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[cc_part]
int Branches::nphe(int i) {
  try {
    return _nphe[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[cc_part]
float Branches::cc_t(int i) {
  try {
    return _cc_t[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[cc_part]
float Branches::cc_r(int i) {
  try {
    return _cc_r[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[cc_part]
float Branches::cc_c2(int i) {
  try {
    return _cc_c2[i];
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}  //[cc_part]

std::vector<int> Branches::id() {
  // [gpart]
  std::vector<int> v;
  for (int i = 0; i < _gpart; i++) v.emplace_back(_id[i]);
  return v;
}
std::vector<int> Branches::stat() {
  // [gpart]
  std::vector<int> v;
  for (int i = 0; i < _gpart; i++) v.emplace_back(_stat[i]);
  return v;
}
std::vector<int> Branches::dc() {
  // [gpart]
  std::vector<int> v;
  for (int i = 0; i < _gpart; i++) v.emplace_back(_dc[i]);
  return v;
}
std::vector<int> Branches::cc() {
  // [gpart]
  std::vector<int> v;
  for (int i = 0; i < _gpart; i++) v.emplace_back(_cc[i]);
  return v;
}
std::vector<int> Branches::sc() {
  // [gpart]
  std::vector<int> v;
  for (int i = 0; i < _gpart; i++) v.emplace_back(_sc[i]);
  return v;
}
std::vector<int> Branches::ec() {
  // [gpart]
  std::vector<int> v;
  for (int i = 0; i < _gpart; i++) v.emplace_back(_ec[i]);
  return v;
}
std::vector<int> Branches::lec() {
  // [gpart]
  std::vector<int> v;
  for (int i = 0; i < _gpart; i++) v.emplace_back(_lec[i]);
  return v;
}
std::vector<int> Branches::ccst() {
  // [gpart]
  std::vector<int> v;
  for (int i = 0; i < _gpart; i++) v.emplace_back(_ccst[i]);
  return v;
}
std::vector<float> Branches::p() {
  // [gpart]
  std::vector<float> v;
  for (int i = 0; i < _gpart; i++) v.emplace_back(_p[i]);
  return v;
}
std::vector<float> Branches::m() {
  // [gpart]
  std::vector<float> v;
  for (int i = 0; i < _gpart; i++) v.emplace_back(_m[i]);
  return v;
}
std::vector<int> Branches::q() {
  // [gpart]
  std::vector<int> v;
  for (int i = 0; i < _gpart; i++) v.emplace_back(_q[i]);
  return v;
}
std::vector<float> Branches::b() {
  // [gpart]
  std::vector<float> v;
  for (int i = 0; i < _gpart; i++) v.emplace_back(_b[i]);
  return v;
}
std::vector<float> Branches::cx() {
  // [gpart]
  std::vector<float> v;
  for (int i = 0; i < _gpart; i++) v.emplace_back(_cx[i]);
  return v;
}
std::vector<float> Branches::cy() {
  // [gpart]
  std::vector<float> v;
  for (int i = 0; i < _gpart; i++) v.emplace_back(_cy[i]);
  return v;
}
std::vector<float> Branches::cz() {
  // [gpart]
  std::vector<float> v;
  for (int i = 0; i < _gpart; i++) v.emplace_back(_cz[i]);
  return v;
}
std::vector<float> Branches::vx() {
  // [gpart]
  std::vector<float> v;
  for (int i = 0; i < _gpart; i++) v.emplace_back(_vx[i]);
  return v;
}
std::vector<float> Branches::vy() {
  // [gpart]
  std::vector<float> v;
  for (int i = 0; i < _gpart; i++) v.emplace_back(_vy[i]);
  return v;
}
std::vector<float> Branches::vz() {
  // [gpart]
  std::vector<float> v;
  for (int i = 0; i < _gpart; i++) v.emplace_back(_vz[i]);
  return v;
}

/*
int Branches::lec_sect(i) { return _; }    //[lac_part]
int Branches::lec_hit(i) { return _; }     //[lac_part]
int Branches::lec_stat(i) { return _; }    //[lac_part]
float Branches::lec_etot(i) { return _; }  //[lac_part]
float Branches::lec_ein(i) { return _; }   //[lac_part]
float Branches::lec_t(i) { return _; }     //[lac_part]
float Branches::lec_r(i) { return _; }     //[lac_part]
float Branches::lec_x(i) { return _; }     //[lac_part]
float Branches::lec_y(i) { return _; }     //[lac_part]
float Branches::lec_z(i) { return _; }     //[lac_part]
float Branches::lec_c2(i) { return _; }    //[lac_part]

int Branches::st_status(i) { return _; }  //[st_part]
float Branches::st_time(i) { return _; }  //[st_part]
float Branches::st_rtrk(i) { return _; }  //[st_part]
*/
