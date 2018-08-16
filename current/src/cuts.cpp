/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include <iostream>
#include "cuts.hpp"

Cuts::Cuts() {}
Cuts::~Cuts() {}

void Cuts::Set_num_phe(int set) { num_phe = set; }
void Cuts::Set_charge(int set) { charge = set; }
void Cuts::Set_electron_id(int set) { electron_id = set; }
void Cuts::Set_gpart(int set) { gpart = set; }
void Cuts::Set_BeamPosition(double x, double y, double z) {
  _vx = x;
  _vy = y;
  _vz = z;
}
void Cuts::Set_elec_fid(double theta, double phi, int sec) {
  _theta = theta;
  _phi = phi;
  _sec = sec + 1;

  switch (_sec) {
    case 1:
      _phi_cent = _phi - 30;
      break;
    case 2:
      _phi_cent = _phi - 90;
      break;
    case 3:
      _phi_cent = _phi - 150;
      break;
    case 4:
      _phi_cent = _phi + 150;
      break;
    case 5:
      _phi_cent = _phi + 90;
      break;
    case 6:
      _phi_cent = _phi + 30;
      break;
  }
}
void Cuts::Set_p(double set) { electron_p = set; }
void Cuts::Set_Sf(double set) { samp_frac = set; }

void Cuts::Set_ec_cut(bool set) { ec_cut = set; }
void Cuts::Set_cc_cut(bool set) { cc_cut = set; }
void Cuts::Set_stat_cut(bool set) { stat_cut = set; }
void Cuts::Set_sc_cut(bool set) { sc_cut = set; }
void Cuts::Set_dc_cut(bool set) { dc_cut = set; }
void Cuts::Set_dc_stat_cut(bool set) { dc_stat_cut = set; }

bool Cuts::isElecctron() {
  bool _elec = true;
  //_elec &= (charge == NEGATIVE);
  _elec &= (gpart > 1);
  //_elec &= ec_cut;
  //_elec &= cc_cut;
  //_elec &= stat_cut;
  //_elec &= sc_cut;
  //_elec &= dc_cut;
  //_elec &= dc_stat_cut;

  electron_cut = _elec;
  return _elec;
}

bool Cuts::isStrictElecctron() {
  bool _elec = true;
  _elec &= isElecctron();
  //_elec &= (electron_p > MIN_P_CUT);

  //_elec &= (num_phe > 20);
  //_elec &= sf_cut(samp_frac, electron_p);
  _elec &= elec_fid_cut();

  //_elec &= (_vx > 0.2);
  //_elec &= (_vx < 0.4);
  //_elec &= (abs(_vy) < 0.2);
  //_elec &= (abs(_vz) < 3);

  electron_cut = _elec;
  return _elec;
}

bool Cuts::CheckElectron() { return electron_cut; }

double Cuts::sf_top_fit(double P) {
  // Old fit for all p...
  // double par[3] = {0.363901, -0.00992778, 5.84749e-06};
  double par[3] = {0.335, 0.003358, -2.551e-6};
  double x[1] = {P};
  return func::ec_fit_func(x, par);
}
double Cuts::sf_bot_fit(double P) {
  // Old fit for all p...
  // double par[3] = {0.103964, 0.0524214, -3.64355e-05};
  double par[3] = {0.1542, 0.02947, -2.411e-5};
  double x[1] = {P};
  return func::ec_fit_func(x, par);
}
bool Cuts::sf_cut(double sf, double P) { return ((sf > sf_bot_fit(P)) && (sf < sf_top_fit(P))); }

double Cuts::dt_P_bot_fit(double P) {
  double par[2] = {-1.509, 0.4172};
  double x[1] = {P};
  return func::dt_fit(x, par);
}
double Cuts::dt_P_top_fit(double P) {
  double par[2] = {1.307, -0.3473};
  double x[1] = {P};
  return func::dt_fit(x, par);
}
bool Cuts::dt_P_cut(double dt, double P) { return ((dt > dt_P_bot_fit(P)) && (dt < dt_P_top_fit(P))); }

double Cuts::dt_Pip_bot_fit(double P) {
  double par[2] = {-0.9285, -0.04094};
  double x[1] = {P};
  return func::dt_fit(x, par);
}
double Cuts::dt_Pip_top_fit(double P) {
  double par[2] = {0.9845, -0.05473};
  double x[1] = {P};
  return func::dt_fit(x, par);
}
bool Cuts::dt_Pip_cut(double dt, double P) { return (dt > dt_Pip_bot_fit(P)) && (dt < dt_Pip_top_fit(P)); }

bool Cuts::elec_fid_cut() {
  double c[4] = {0.0548311203, 0.0327878012, 17.0683287374};
  double y = c[0] * _phi_cent * _phi_cent + c[1] * _phi_cent + c[2];

  return _theta >= y;
}
