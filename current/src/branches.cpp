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

Branches::Branches(TChain* tree, bool MC) {
  _MC = MC;
  myTree = tree;
  init();
}

Branches::Branches(const Branches& b) {
  myTree = b.myTree;
  _MC = b._MC;
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
  if (_MC) {
    myTree->SetBranchAddress("nprt", &_nprt);
    myTree->SetBranchAddress("pidpart", _pidpart);
    myTree->SetBranchAddress("xpart", _xpart);
    myTree->SetBranchAddress("ypart", _ypart);
    myTree->SetBranchAddress("zpart", _zpart);
    myTree->SetBranchAddress("epart", _epart);
    myTree->SetBranchAddress("pxpart", _pxpart);
    myTree->SetBranchAddress("pypart", _pypart);
    myTree->SetBranchAddress("pzpart", _pzpart);
    myTree->SetBranchAddress("qpart", _qpart);
  }
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

int Branches::nprt() { return _nprt; }

int Branches::id(int i) {
  if (i < _gpart) {
    return _id[i];
  } else {
    return int(NULL);
  }
}  // [gpart]
int Branches::stat(int i) {
  if (i < _gpart) {
    return _stat[i];
  } else {
    return int(NULL);
  }
}  // [gpart]
int Branches::dc(int i) {
  if (i < _gpart) {
    return _dc[i];
  } else {
    return int(NULL);
  }
}  // [gpart]
int Branches::cc(int i) {
  if (i < _gpart) {
    return _cc[i];
  } else {
    return int(NULL);
  }
}  // [gpart]
int Branches::sc(int i) {
  if (i < _gpart) {
    return _sc[i];
  } else {
    return int(NULL);
  }
}  // [gpart]
int Branches::ec(int i) {
  if (i < _gpart) {
    return _ec[i];
  } else {
    return int(NULL);
  }
}  // [gpart]
int Branches::lec(int i) {
  if (i < _gpart) {
    return _lec[i];
  } else {
    return int(NULL);
  }
}  // [gpart]
int Branches::ccst(int i) {
  if (i < _gpart) {
    return _ccst[i];
  } else {
    return int(NULL);
  }
}  // [gpart]
float Branches::p(int i) {
  if (i < _gpart) {
    return _p[i];
  } else {
    return std::nanf("NULL");
  }
}  // [gpart]
float Branches::px(int i) {
  if (i < _gpart) {
    return _p[i] * _cx[i];
  } else {
    return std::nanf("NULL");
  }
}  // [gpart]
float Branches::py(int i) {
  if (i < _gpart) {
    return _p[i] * _cy[i];
  } else {
    return std::nanf("NULL");
  }
}  // [gpart]
float Branches::pz(int i) {
  if (i < _gpart) {
    return _p[i] * _cz[i];
  } else {
    return std::nanf("NULL");
  }
}  // [gpart]
float Branches::m(int i) {
  if (i < _gpart) {
    return _m[i];
  } else {
    return std::nanf("NULL");
  }
}  // [gpart]
int Branches::q(int i) {
  if (i < _gpart) {
    return _q[i];
  } else {
    return int(NULL);
  }
}  // [gpart]
float Branches::b(int i) {
  if (i < _gpart) {
    return _b[i];
  } else {
    return std::nanf("NULL");
  }
}  // [gpart]
float Branches::cx(int i) {
  if (i < _gpart) {
    return _cx[i];
  } else {
    return std::nanf("NULL");
  }
}  // [gpart]
float Branches::cy(int i) {
  if (i < _gpart) {
    return _cy[i];
  } else {
    return std::nanf("NULL");
  }
}  // [gpart]
float Branches::cz(int i) {
  if (i < _gpart) {
    return _cz[i];
  } else {
    return std::nanf("NULL");
  }
}  // [gpart]
float Branches::vx(int i) {
  if (i < _gpart) {
    return _vx[i];
  } else {
    return std::nanf("NULL");
  }
}  // [gpart]
float Branches::vy(int i) {
  if (i < _gpart) {
    return _vy[i];
  } else {
    return std::nanf("NULL");
  }
}  // [gpart]
float Branches::vz(int i) {
  if (i < _gpart) {
    return _vz[i];
  } else {
    return std::nanf("NULL");
  }
}  // [gpart]

int Branches::dc_sect(int i) {
  if (i < _dc_part) {
    return _dc_sect[_dc[i] - 1];
  } else {
    return int(NULL);
  }
}  //[dc_part]
int Branches::dc_trk(int i) {
  if (i < _dc_part) {
    return _dc_trk[_dc[i] - 1];
  } else {
    return int(NULL);
  }
}  //[dc_part]
int Branches::dc_stat(int i) {
  if (i < _dc_part) {
    return _dc_stat[_dc[i] - 1];
  } else {
    return int(NULL);
  }
}  //[dc_part]
float Branches::dc_vx(int i) {
  if (i < _dc_part) {
    return _dc_vx[_dc[i] - 1];
  } else {
    return std::nanf("NULL");
  }
}  //[dc_part]
float Branches::dc_vy(int i) {
  if (i < _dc_part) {
    return _dc_vy[_dc[i] - 1];
  } else {
    return std::nanf("NULL");
  }
}  //[dc_part]
float Branches::dc_vz(int i) {
  if (i < _dc_part) {
    return _dc_vz[_dc[i] - 1];
  } else {
    return std::nanf("NULL");
  }
}  //[dc_part]
float Branches::dc_vr(int i) {
  if (i < _dc_part) {
    return _dc_vr[_dc[i] - 1];
  } else {
    return std::nanf("NULL");
  }
}  //[dc_part]
float Branches::dc_xsc(int i) {
  if (i < _dc_part) {
    return _dc_xsc[_dc[i] - 1];
  } else {
    return std::nanf("NULL");
  }
}  //[dc_part]
float Branches::dc_ysc(int i) {
  if (i < _dc_part) {
    return _dc_ysc[_dc[i] - 1];
  } else {
    return std::nanf("NULL");
  }
}  //[dc_part]
float Branches::dc_zsc(int i) {
  if (i < _dc_part) {
    return _dc_zsc[_dc[i] - 1];
  } else {
    return std::nanf("NULL");
  }
}  //[dc_part]
float Branches::dc_cxsc(int i) {
  if (i < _dc_part) {
    return _dc_cxsc[_dc[i] - 1];
  } else {
    return std::nanf("NULL");
  }
}  //[dc_part]
float Branches::dc_cysc(int i) {
  if (i < _dc_part) {
    return _dc_cysc[_dc[i] - 1];
  } else {
    return std::nanf("NULL");
  }
}  //[dc_part]
float Branches::dc_czsc(int i) {
  if (i < _dc_part) {
    return _dc_czsc[_dc[i] - 1];
  } else {
    return std::nanf("NULL");
  }
}  //[dc_part]
float Branches::dc_c2(int i) {
  if (i < _dc_part) {
    return _dc_c2[_dc[i] - 1];
  } else {
    return std::nanf("NULL");
  }
}  //[dc_part]

int Branches::ec_stat(int i) {
  if (i < _ec_part) {
    return _ec_stat[_ec[i] - 1];
  } else {
    return int(NULL);
  }
}  //[ec_part]
int Branches::ec_sect(int i) {
  if (i < _ec_part) {
    return _ec_sect[_ec[i] - 1];
  } else {
    return int(NULL);
  }
}  //[ec_part]
int Branches::ec_whol(int i) {
  if (i < _ec_part) {
    return _ec_whol[_ec[i] - 1];
  } else {
    return int(NULL);
  }
}  //[ec_part]
int Branches::ec_inst(int i) {
  if (i < _ec_part) {
    return _ec_inst[_ec[i] - 1];
  } else {
    return int(NULL);
  }
}  //[ec_part]
int Branches::ec_oust(int i) {
  if (i < _ec_part) {
    return _ec_oust[_ec[i] - 1];
  } else {
    return int(NULL);
  }
}  //[ec_part]
float Branches::etot(int i) {
  if (i < _ec_part) {
    return _etot[_ec[i] - 1];
  } else {
    return std::nanf("NULL");
  }
}  //[ec_part]
float Branches::ec_ei(int i) {
  if (i < _ec_part) {
    return _ec_ei[_ec[i] - 1];
  } else {
    return std::nanf("NULL");
  }
}  //[ec_part]
float Branches::ec_eo(int i) {
  if (i < _ec_part) {
    return _ec_eo[_ec[i] - 1];
  } else {
    return std::nanf("NULL");
  }
}  //[ec_part]
float Branches::ec_t(int i) {
  if (i < _ec_part) {
    return _ec_t[_ec[i] - 1];
  } else {
    return std::nanf("NULL");
  }
}  //[ec_part]
float Branches::ec_r(int i) {
  if (i < _ec_part) {
    return _ec_r[_ec[i] - 1];
  } else {
    return std::nanf("NULL");
  }
}  //[ec_part]
float Branches::ech_x(int i) {
  if (i < _ec_part) {
    return _ech_x[_ec[i] - 1];
  } else {
    return std::nanf("NULL");
  }
}  //[ec_part]
float Branches::ech_y(int i) {
  if (i < _ec_part) {
    return _ech_y[_ec[i] - 1];
  } else {
    return std::nanf("NULL");
  }
}  //[ec_part]
float Branches::ech_z(int i) {
  if (i < _ec_part) {
    return _ech_z[_ec[i] - 1];
  } else {
    return std::nanf("NULL");
  }
}  //[ec_part]
float Branches::ec_m2(int i) {
  if (i < _ec_part) {
    return _ec_m2[_ec[i] - 1];
  } else {
    return std::nanf("NULL");
  }
}  //[ec_part]
float Branches::ec_m3(int i) {
  if (i < _ec_part) {
    return _ec_m3[_ec[i] - 1];
  } else {
    return std::nanf("NULL");
  }
}  //[ec_part]
float Branches::ec_m4(int i) {
  if (i < _ec_part) {
    return _ec_m4[_ec[i] - 1];
  } else {
    return std::nanf("NULL");
  }
}  //[ec_part]
float Branches::ec_c2(int i) {
  if (i < _ec_part) {
    return _ec_c2[_ec[i] - 1];
  } else {
    return std::nanf("NULL");
  }
}  //[ec_part]

int Branches::sc_sect(int i) {
  if (i < _sc_part) {
    return _sc_sect[_sc[i] - 1];
  } else {
    return int(NULL);
  }
}  //[sc_part]
int Branches::sc_hit(int i) {
  if (i < _sc_part) {
    return _sc_hit[_sc[i] - 1];
  } else {
    return int(NULL);
  }
}  //[sc_part]
int Branches::sc_pd(int i) {
  if (i < _sc_part) {
    return _sc_pd[_sc[i] - 1];
  } else {
    return int(NULL);
  }
}  //[sc_part]
int Branches::sc_stat(int i) {
  if (i < _sc_part) {
    return _sc_stat[_sc[i] - 1];
  } else {
    return int(NULL);
  }
}  //[sc_part]
float Branches::edep(int i) {
  if (i < _sc_part) {
    return _edep[_sc[i] - 1];
  } else {
    return std::nanf("NULL");
  }
}  //[sc_part]
float Branches::sc_t(int i) {
  if (i < _sc_part) {
    return _sc_t[_sc[i] - 1];
  } else {
    return std::nanf("NULL");
  }
}  //[sc_part]
float Branches::sc_r(int i) {
  if (i < _sc_part) {
    return _sc_r[_sc[i] - 1];
  } else {
    return std::nanf("NULL");
  }
}  //[sc_part]
float Branches::sc_c2(int i) {
  if (i < _sc_part) {
    return _sc_c2[_sc[i] - 1];
  } else {
    return std::nanf("NULL");
  }
}  //[sc_part]

int Branches::cc_sect(int i) {
  if (i < _cc_part) {
    return _cc_sect[_cc[i] - 1];
  } else {
    return int(NULL);
  }
}  //[cc_part]
int Branches::cc_hit(int i) {
  if (i < _cc_part) {
    return _cc_hit[_cc[i] - 1];
  } else {
    return int(NULL);
  }
}  //[cc_part]
int Branches::cc_segm(int i) {
  if (i < _cc_part) {
    return _cc_segm[_cc[i] - 1];
  } else {
    return int(NULL);
  }
}  //[cc_part]
int Branches::nphe(int i) {
  if (i < _cc_part) {
    return _nphe[_cc[i] - 1];
  } else {
    return int(NULL);
  }
}  //[cc_part]
float Branches::cc_t(int i) {
  if (i < _cc_part) {
    return _cc_t[_cc[i] - 1];
  } else {
    return std::nanf("NULL");
  }
}  //[cc_part]
float Branches::cc_r(int i) {
  if (i < _cc_part) {
    return _cc_r[_cc[i] - 1];
  } else {
    return std::nanf("NULL");
  }
}  //[cc_part]
float Branches::cc_c2(int i) {
  if (i < _cc_part) {
    return _cc_c2[_cc[i] - 1];
  } else {
    return std::nanf("NULL");
  }
}  //[cc_part]

int Branches::pidpart(int i) {
  if (i < _nprt) {
    return _pidpart[i];
  } else {
    return (int)NULL;
  }
}  //[nprt]

float Branches::xpart(int i) {
  if (i < _nprt) {
    return _xpart[i];
  } else {
    return (float)NULL;
  }
}  //[nprt]

float Branches::ypart(int i) {
  if (i < _nprt) {
    return _ypart[i];
  } else {
    return (float)NULL;
  }
}  //[nprt]

float Branches::zpart(int i) {
  if (i < _nprt) {
    return _zpart[i];
  } else {
    return (float)NULL;
  }
}  //[nprt]

float Branches::epart(int i) {
  if (i < _nprt) {
    return _epart[i];
  } else {
    return (float)NULL;
  }
}  //[nprt]

float Branches::pxpart(int i) {
  if (i < _nprt) {
    return _pxpart[i];
  } else {
    return (float)NULL;
  }
}  //[nprt]

float Branches::pypart(int i) {
  if (i < _nprt) {
    return _pypart[i];
  } else {
    return (float)NULL;
  }
}  //[nprt]

float Branches::pzpart(int i) {
  if (i < _nprt) {
    return _pzpart[i];
  } else {
    return (float)NULL;
  }
}  //[nprt]

float Branches::qpart(int i) {
  if (i < _nprt) {
    return _qpart[i];
  } else {
    return (float)NULL;
  }
}  //[nprt]

std::vector<int> Branches::id() {
  // [gpart]
  std::vector<int> v(_gpart);
  for (int i = 0; i < _gpart; i++) v[i] = _id[i];
  return v;
}
std::vector<int> Branches::stat() {
  // [gpart]
  std::vector<int> v(_gpart);
  for (int i = 0; i < _gpart; i++) v[i] = _stat[i];
  return v;
}
std::vector<int> Branches::dc() {
  // [gpart]
  std::vector<int> v(_gpart);
  for (int i = 0; i < _gpart; i++) v[i] = _dc[i];
  return v;
}
std::vector<int> Branches::cc() {
  // [gpart]
  std::vector<int> v(_gpart);
  for (int i = 0; i < _gpart; i++) v[i] = _cc[i];
  return v;
}
std::vector<int> Branches::sc() {
  // [gpart]
  std::vector<int> v(_gpart);
  for (int i = 0; i < _gpart; i++) v[i] = _sc[i];
  return v;
}
std::vector<int> Branches::ec() {
  // [gpart]
  std::vector<int> v(_gpart);
  for (int i = 0; i < _gpart; i++) v[i] = _ec[i];
  return v;
}
std::vector<int> Branches::lec() {
  // [gpart]
  std::vector<int> v(_gpart);
  for (int i = 0; i < _gpart; i++) v[i] = _lec[i];
  return v;
}
std::vector<int> Branches::ccst() {
  // [gpart]
  std::vector<int> v(_gpart);
  for (int i = 0; i < _gpart; i++) v[i] = _ccst[i];
  return v;
}
std::vector<float> Branches::p() {
  // [gpart]
  std::vector<float> v(_gpart);
  for (int i = 0; i < _gpart; i++) v[i] = _p[i];
  return v;
}
std::vector<float> Branches::px() {
  // [gpart]
  std::vector<float> v(_gpart);
  for (int i = 0; i < _gpart; i++) v[i] = _p[i] * _cx[i];
  return v;
}
std::vector<float> Branches::py() {
  // [gpart]
  std::vector<float> v(_gpart);
  for (int i = 0; i < _gpart; i++) v[i] = _p[i] * _cy[i];
  return v;
}
std::vector<float> Branches::pz() {
  // [gpart]
  std::vector<float> v(_gpart);
  for (int i = 0; i < _gpart; i++) v[i] = _p[i] * _cz[i];
  return v;
}
std::vector<float> Branches::m() {
  // [gpart]
  std::vector<float> v(_gpart);
  for (int i = 0; i < _gpart; i++) v[i] = _m[i];
  return v;
}
std::vector<int> Branches::q() {
  // [gpart]
  std::vector<int> v(_gpart);
  for (int i = 0; i < _gpart; i++) v[i] = _q[i];
  return v;
}
std::vector<float> Branches::b() {
  // [gpart]
  std::vector<float> v(_gpart);
  for (int i = 0; i < _gpart; i++) v[i] = _b[i];
  return v;
}
std::vector<float> Branches::cx() {
  // [gpart]
  std::vector<float> v(_gpart);
  for (int i = 0; i < _gpart; i++) v[i] = _cx[i];
  return v;
}
std::vector<float> Branches::cy() {
  // [gpart]
  std::vector<float> v(_gpart);
  for (int i = 0; i < _gpart; i++) v[i] = _cy[i];
  return v;
}
std::vector<float> Branches::cz() {
  // [gpart]
  std::vector<float> v(_gpart);
  for (int i = 0; i < _gpart; i++) v[i] = _cz[i];
  return v;
}
std::vector<float> Branches::vx() {
  // [gpart]
  std::vector<float> v(_gpart);
  for (int i = 0; i < _gpart; i++) v[i] = _vx[i];
  return v;
}
std::vector<float> Branches::vy() {
  // [gpart]
  std::vector<float> v(_gpart);
  for (int i = 0; i < _gpart; i++) v[i] = _vy[i];
  return v;
}
std::vector<float> Branches::vz() {
  // [gpart]
  std::vector<float> v(_gpart);
  for (int i = 0; i < _gpart; i++) v[i] = _vz[i];
  return v;
}

std::vector<int> Branches::dc_sect() {
  // [dcpart]
  std::vector<int> v(_dc_part);
  for (int i = 0; i < _dc_part; i++) v[i] = _dc_sect[i];
  return v;
}
std::vector<int> Branches::dc_trk() {
  // [dcpart]
  std::vector<int> v(_dc_part);
  for (int i = 0; i < _dc_part; i++) v[i] = _dc_trk[i];
  return v;
}
std::vector<int> Branches::dc_stat() {
  // [dcpart]
  std::vector<int> v(_dc_part);
  for (int i = 0; i < _dc_part; i++) v[i] = _dc_stat[i];
  return v;
}
std::vector<float> Branches::dc_vx() {
  // [dcpart]
  std::vector<float> v(_dc_part);
  for (int i = 0; i < _dc_part; i++) v[i] = _dc_vx[i];
  return v;
}
std::vector<float> Branches::dc_vy() {
  // [dcpart]
  std::vector<float> v(_dc_part);
  for (int i = 0; i < _dc_part; i++) v[i] = _dc_vy[i];
  return v;
}
std::vector<float> Branches::dc_vz() {
  // [dcpart]
  std::vector<float> v(_dc_part);
  for (int i = 0; i < _dc_part; i++) v[i] = _dc_vz[i];
  return v;
}
std::vector<float> Branches::dc_vr() {
  // [dcpart]
  std::vector<float> v(_dc_part);
  for (int i = 0; i < _dc_part; i++) v[i] = _dc_vr[i];
  return v;
}
std::vector<float> Branches::dc_xsc() {
  // [dcpart]
  std::vector<float> v(_dc_part);
  for (int i = 0; i < _dc_part; i++) v[i] = _dc_xsc[i];
  return v;
}
std::vector<float> Branches::dc_ysc() {
  // [dcpart]
  std::vector<float> v(_dc_part);
  for (int i = 0; i < _dc_part; i++) v[i] = _dc_ysc[i];
  return v;
}
std::vector<float> Branches::dc_zsc() {
  // [dcpart]
  std::vector<float> v(_dc_part);
  for (int i = 0; i < _dc_part; i++) v[i] = _dc_zsc[i];
  return v;
}
std::vector<float> Branches::dc_cxsc() {
  // [dcpart]
  std::vector<float> v(_dc_part);
  for (int i = 0; i < _dc_part; i++) v[i] = _dc_cxsc[i];
  return v;
}
std::vector<float> Branches::dc_cysc() {
  // [dcpart]
  std::vector<float> v(_dc_part);
  for (int i = 0; i < _dc_part; i++) v[i] = _dc_cysc[i];
  return v;
}
std::vector<float> Branches::dc_czsc() {
  // [dcpart]
  std::vector<float> v(_dc_part);
  for (int i = 0; i < _dc_part; i++) v[i] = _dc_czsc[i];
  return v;
}
std::vector<float> Branches::dc_c2() {
  // [dcpart]
  std::vector<float> v(_dc_part);
  for (int i = 0; i < _dc_part; i++) v[i] = _dc_c2[i];
  return v;
}

std::vector<int> Branches::ec_stat() {
  // [ec_part]
  std::vector<int> v(_ec_part);
  for (int i = 0; i < _ec_part; i++) v[i] = _ec_stat[i];
  return v;
}
std::vector<int> Branches::ec_sect() {
  // [ec_part]
  std::vector<int> v(_ec_part);
  for (int i = 0; i < _ec_part; i++) v[i] = _ec_sect[i];
  return v;
}
std::vector<int> Branches::ec_whol() {
  // [ec_part]
  std::vector<int> v(_ec_part);
  for (int i = 0; i < _ec_part; i++) v[i] = _ec_whol[i];
  return v;
}
std::vector<int> Branches::ec_inst() {
  // [ec_part]
  std::vector<int> v(_ec_part);
  for (int i = 0; i < _ec_part; i++) v[i] = _ec_inst[i];
  return v;
}
std::vector<int> Branches::ec_oust() {
  // [ec_part]
  std::vector<int> v(_ec_part);
  for (int i = 0; i < _ec_part; i++) v[i] = _ec_oust[i];
  return v;
}
std::vector<float> Branches::etot() {
  // [ec_part]
  std::vector<float> v(_ec_part);
  for (int i = 0; i < _ec_part; i++) v[i] = _etot[i];
  return v;
}
std::vector<float> Branches::ec_ei() {
  // [ec_part]
  std::vector<float> v(_ec_part);
  for (int i = 0; i < _ec_part; i++) v[i] = _ec_ei[i];
  return v;
}
std::vector<float> Branches::ec_eo() {
  // [ec_part]
  std::vector<float> v(_ec_part);
  for (int i = 0; i < _ec_part; i++) v[i] = _ec_eo[i];
  return v;
}
std::vector<float> Branches::ec_t() {
  // [ec_part]
  std::vector<float> v(_ec_part);
  for (int i = 0; i < _ec_part; i++) v[i] = _ec_t[i];
  return v;
}
std::vector<float> Branches::ec_r() {
  // [ec_part]
  std::vector<float> v(_ec_part);
  for (int i = 0; i < _ec_part; i++) v[i] = _ec_r[i];
  return v;
}
std::vector<float> Branches::ech_x() {
  // [ec_part]
  std::vector<float> v(_ec_part);
  for (int i = 0; i < _ec_part; i++) v[i] = _ech_x[i];
  return v;
}
std::vector<float> Branches::ech_y() {
  // [ec_part]
  std::vector<float> v(_ec_part);
  for (int i = 0; i < _ec_part; i++) v[i] = _ech_y[i];
  return v;
}
std::vector<float> Branches::ech_z() {
  // [ec_part]
  std::vector<float> v(_ec_part);
  for (int i = 0; i < _ec_part; i++) v[i] = _ech_z[i];
  return v;
}
std::vector<float> Branches::ec_m2() {
  // [ec_part]
  std::vector<float> v(_ec_part);
  for (int i = 0; i < _ec_part; i++) v[i] = _ec_m2[i];
  return v;
}
std::vector<float> Branches::ec_m3() {
  // [ec_part]
  std::vector<float> v(_ec_part);
  for (int i = 0; i < _ec_part; i++) v[i] = _ec_m3[i];
  return v;
}
std::vector<float> Branches::ec_m4() {
  // [ec_part]
  std::vector<float> v(_ec_part);
  for (int i = 0; i < _ec_part; i++) v[i] = _ec_m4[i];
  return v;
}
std::vector<float> Branches::ec_c2() {
  // [ec_part]
  std::vector<float> v(_ec_part);
  for (int i = 0; i < _ec_part; i++) v[i] = _ec_c2[i];
  return v;
}

std::vector<int> Branches::sc_sect() {
  // [sc_part]
  std::vector<int> v(_sc_part);
  for (int i = 0; i < _sc_part; i++) v[i] = _sc_sect[i];
  return v;
}
std::vector<int> Branches::sc_hit() {
  // [sc_part]
  std::vector<int> v(_sc_part);
  for (int i = 0; i < _sc_part; i++) v[i] = _sc_hit[i];
  return v;
}
std::vector<int> Branches::sc_pd() {
  // [sc_part]
  std::vector<int> v(_sc_part);
  for (int i = 0; i < _sc_part; i++) v[i] = _sc_pd[i];
  return v;
}
std::vector<int> Branches::sc_stat() {
  // [sc_part]
  std::vector<int> v(_sc_part);
  for (int i = 0; i < _sc_part; i++) v[i] = _sc_stat[i];
  return v;
}
std::vector<float> Branches::edep() {
  // [sc_part]
  std::vector<float> v(_sc_part);
  for (int i = 0; i < _sc_part; i++) v[i] = _edep[i];
  return v;
}
std::vector<float> Branches::sc_t() {
  // [sc_part]
  std::vector<float> v(_sc_part);
  for (int i = 0; i < _sc_part; i++) v[i] = _sc_t[i];
  return v;
}
std::vector<float> Branches::sc_r() {
  // [sc_part]
  std::vector<float> v(_sc_part);
  for (int i = 0; i < _sc_part; i++) v[i] = _sc_r[i];
  return v;
}
std::vector<float> Branches::sc_c2() {
  // [sc_part]
  std::vector<float> v(_sc_part);
  for (int i = 0; i < _sc_part; i++) v[i] = _sc_c2[i];
  return v;
}

std::vector<int> Branches::cc_sect() {
  // [cc_part]
  std::vector<int> v(_cc_part);
  for (int i = 0; i < _cc_part; i++) v[i] = _cc_sect[i];
  return v;
}
std::vector<int> Branches::cc_hit() {
  // [cc_part]
  std::vector<int> v(_cc_part);
  for (int i = 0; i < _cc_part; i++) v[i] = _cc_hit[i];
  return v;
}
std::vector<int> Branches::cc_segm() {
  // [cc_part]
  std::vector<int> v(_cc_part);
  for (int i = 0; i < _cc_part; i++) v[i] = _cc_segm[i];
  return v;
}
std::vector<int> Branches::nphe() {
  // [cc_part]
  std::vector<int> v(_cc_part);
  for (int i = 0; i < _cc_part; i++) v[i] = _nphe[i];
  return v;
}
std::vector<float> Branches::cc_t() {
  // [cc_part]
  std::vector<float> v(_cc_part);
  for (int i = 0; i < _cc_part; i++) v[i] = _cc_t[i];
  return v;
}
std::vector<float> Branches::cc_r() {
  // [cc_part]
  std::vector<float> v(_cc_part);
  for (int i = 0; i < _cc_part; i++) v[i] = _cc_r[i];
  return v;
}
std::vector<float> Branches::cc_c2() {
  // [cc_part]
  std::vector<float> v(_cc_part);
  for (int i = 0; i < _cc_part; i++) v[i] = _cc_c2[i];
  return v;
}

std::vector<int> Branches::pidpart() {
  //[nprt]
  std::vector<int> v(_nprt);
  for (int i = 0; i < _nprt; i++) v[i] = _pidpart[i];
  return v;
}

std::vector<float> Branches::xpart() {
  //[nprt]
  std::vector<float> v(_nprt);
  for (int i = 0; i < _nprt; i++) v[i] = _xpart[i];
  return v;
}

std::vector<float> Branches::ypart() {
  //[nprt]
  std::vector<float> v(_nprt);
  for (int i = 0; i < _nprt; i++) v[i] = _ypart[i];
  return v;
}

std::vector<float> Branches::zpart() {
  //[nprt]
  std::vector<float> v(_nprt);
  for (int i = 0; i < _nprt; i++) v[i] = _zpart[i];
  return v;
}

std::vector<float> Branches::epart() {
  //[nprt]
  std::vector<float> v(_nprt);
  for (int i = 0; i < _nprt; i++) v[i] = _epart[i];
  return v;
}

std::vector<float> Branches::pxpart() {
  //[nprt]
  std::vector<float> v(_nprt);
  for (int i = 0; i < _nprt; i++) v[i] = _pxpart[i];
  return v;
}

std::vector<float> Branches::pypart() {
  //[nprt]
  std::vector<float> v(_nprt);
  for (int i = 0; i < _nprt; i++) v[i] = _pypart[i];
  return v;
}

std::vector<float> Branches::pzpart() {
  //[nprt]
  std::vector<float> v(_nprt);
  for (int i = 0; i < _nprt; i++) v[i] = _pzpart[i];
  return v;
}

std::vector<float> Branches::qpart() {
  //[nprt]
  std::vector<float> v(_nprt);
  for (int i = 0; i < _nprt; i++) v[i] = _qpart[i];
  return v;
}
