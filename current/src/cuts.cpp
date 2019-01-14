/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include <iostream>
#include "cuts.hpp"
Cuts::Cuts() {}
Cuts::Cuts(Branches* data) { _data = data; }
Cuts::~Cuts() {}

void Cuts::Set_elec_fid() {
  _theta = physics::theta_calc(_data->cz(0));
  _phi = physics::phi_calc(_data->cx(0), _data->cy(0));
  _sector = physics::get_sector(_phi) + 1;

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

bool Cuts::isElecctron() {
  bool _elec = true;
  _elec &= (_data->q(0) == NEGATIVE);
  _elec &= (_data->ec(0) > 0);    // ``` ``` ``` ec
  _elec &= (_data->gpart() > 0);  // Number of good particles is greater than 0
  _elec &= (_data->cc(0) > 0);
  _elec &= (_data->stat(0) > 0);  // First Particle hit stat
  _elec &= (_data->sc(0) > 0);
  _elec &= (_data->dc(0) > 0);
  _elec &= (_data->nphe(0) > 30);
  _elec &= (_data->dc_stat(0) > 0);

  return _elec;
}

bool Cuts::Fid_cut() {
  Set_elec_fid();
  return elec_fid_cut();
}

bool Cuts::Beam_cut() {
  bool _beam = true;
  _beam &= (_data->vx(0) > 0.2);
  _beam &= (_data->vx(0) < 0.4);
  _beam &= (abs(_data->vy(0)) < 0.1);
  _beam &= (abs(_data->vz(0)) < 3);
  return (isElecctron() && _beam);
}

bool Cuts::isStrictElecctron() {
  bool _elec = true;
  _elec &= isElecctron();
  _elec &= Fid_cut();
  _elec &= (electron_p > MIN_P_CUT);

  _elec &= (_data->nphe(0) > 30);
  _elec &= sf_cut(samp_frac, electron_p);

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
