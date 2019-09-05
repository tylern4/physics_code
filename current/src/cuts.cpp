/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "cuts.hpp"
#include <iostream>
//#include "fid_cuts.hpp"

Cuts::Cuts(std::shared_ptr<Branches> data) : _data(data) {}
Cuts::~Cuts() {}

void Cuts::Set_elec_fid() {
  _theta = physics::theta_calc(_data->cz(0));
  _phi = physics::phi_calc(_data->cx(0), _data->cy(0));
  _sec = _data->dc_sect(0);

  switch (_sec) {
    case 1:
      _phi_cent = _phi - 90;
      break;
    case 2:
      _phi_cent = _phi - 30;
      break;
    case 3:
      _phi_cent = _phi + 30;
      break;
    case 4:
      _phi_cent = _phi + 90;
      break;
    case 5:
      _phi_cent = _phi + 150;
      break;
    case 6:
      _phi_cent = _phi - 150;
      break;
  }
}

bool Cuts::isElecctron() {
  bool _elec = true;
  _elec &= (_data->gpart() > 0);  // Number of good particles is greater than 0
  // _elec &= (_data->gpart() < 5);
  if (!_elec) return false;
  _elec &= (_data->q(0) == NEGATIVE);
  _elec &= (_data->ec(0) > 0);
  _elec &= (_data->cc(0) > 0);
  _elec &= (_data->stat(0) > 0);  // First Particle stat
  _elec &= (_data->sc(0) > 0);
  _elec &= (_data->dc(0) > 0);
  _elec &= (_data->dc_stat(0) > 0);

  // Sampling fraction cut
  ////_elec &= sf_cut(_data->etot(0) / _data->p(0), _data->p(0));
  // Cut out low ec inner
  ////_elec &= (_data->ec_ei(0) >= 0.05);
  // Minimum momentum cut
  ////_elec &= (_data->p(0) > MIN_P_CUT);
  // Beam Position cut
  ////_elec &= Beam_cut();
  // Fid Cuts
  ////_elec &= Fid_cut();

  return _elec;
}

bool Cuts::Fid_cut() {
  Set_elec_fid();
  return elec_fid_cut();
}

bool Cuts::Beam_cut() {
  bool _beam = true;

  _beam &= (_data->dc_vx(0) > 0.2 && _data->dc_vx(0) < 0.4);
  _beam &= (_data->dc_vy(0) > -0.1 && _data->dc_vy(0) < 0.16);
  _beam &= (_data->dc_vz(0) > -2.0 && _data->dc_vz(0) < 2.0);

  for (short i = 0; i < _data->gpart(); i++) {
    _beam &= (_data->vx(i) > -1.0 && _data->vx(i) < 1.0);
    _beam &= (_data->vy(i) > -1.5 && _data->vy(i) < 1.5);
    _beam &= (_data->vz(i) > -2.0 && _data->vz(i) < 2.0);
  }

  return _beam;
}

bool Cuts::isStrictElecctron() {
  bool _elec = true;
  _elec &= isElecctron();
  // remove CC hit to both
  ////_elec &= (_data->cc_segm(0) / 1000 - 1 != 0);
  // Cut low number of photo electrons in cc
  ////_elec &= (_data->nphe(0) > 15);

  return _elec;
}

double Cuts::sf_top_fit(double P) {
  double par[3] = {0.368209, 0.000961273, 4.8e-07};
  // double par[3] = {0.3269, 0.000336, 7.731e-7};
  double x[1] = {P};
  return func::ec_fit_func(x, par);
}
double Cuts::sf_bot_fit(double P) {
  double par[3] = {0.162189, 0.0134756, -2e-05};
  // double par[3] = {0.1787, 0.02032, -2.726e-6};
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
  // double par[2] = {-0.9285, -0.04094};
  double par[2] = {-1.3, 0.2616};
  double x[1] = {P};
  return func::dt_fit(x, par);
}
double Cuts::dt_Pip_top_fit(double P) {
  // double par[2] = {0.9845, -0.05473};
  double par[2] = {1.461, -0.4109};
  double x[1] = {P};
  return func::dt_fit(x, par);
}
bool Cuts::dt_Pip_cut(double dt, double P) { return (dt > dt_Pip_bot_fit(P)) && (dt < dt_Pip_top_fit(P)); }

bool Cuts::elec_fid_cut() {
  // double c[3] = {0.0548311203, 0.0327878012, 17.0683287374};
  double c[3] = {0.04, 0.03, 15};
  double y = c[0] * _phi_cent * _phi_cent + c[1] * _phi_cent + c[2];

  return _theta >= y;
}
