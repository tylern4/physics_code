/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "cuts.hpp"
#include "fid_cuts.hpp"

Cuts::Cuts(const std::shared_ptr<Branches>& data) : _data(data) {
  if (_data->gpart() > 0) _dt = std::make_shared<Delta_T>(_data);
}
Cuts::Cuts(const std::shared_ptr<Branches>& data, const std::shared_ptr<Delta_T>& dt) : _data(data), _dt(dt) {}
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

float Cuts::hardon_fid_phi(int part) {
  float phi = physics::phi_calc(_data->cx(part), _data->cy(part));
  short sec = _data->dc_sect(part);

  switch (sec) {
    case 1:
      return phi - 90;
    case 2:
      return phi - 30;
    case 3:
      return phi + 30;
    case 4:
      return phi + 90;
    case 5:
      return phi + 150;
    case 6:
      return phi - 150;
  }
  return NAN;
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
  fid_e(_data->p(0), _data->cz(0), _data->cx(0), _data->cy(0));
  // Sampling fraction cut
  _elec &= sf_cut(_data->etot(0) / _data->p(0), _data->p(0));
  // Cut out low ec inner
  _elec &= (_data->ec_ei(0) >= 0.05);
  // Minimum momentum cut
  //_elec &= (_data->p(0) > MIN_P_CUT);
  // Beam Position cut
  _elec &= Beam_cut();
  // Fid Cuts
  //_elec &= Fid_cut();

  return _elec;
}

float Cuts::hadron_fid_phi_min(float theta, int idx) {
  return -(a0mh[idx] * (1.0 - expf(-a1mh[idx] * (theta - a2mh[idx]))) - a3mh[idx]);
}

float Cuts::hadron_fid_phi_max(float theta, int idx) {
  return (a0xh[idx] * (1.0 - expf(-a1xh[idx] * (theta - a2xh[idx]))) + a3xh[idx]);
}

bool Cuts::Hardon_fid_arjun(int part) {
  float theta = physics::theta_calc(_data->cz(part));
  float phi_c = hardon_fid_phi(part);
  int sector = _data->dc_sect(part);

  if (sector == 0) return false;

  if (phi_c >= hadron_fid_phi_min(theta, sector - 1) && phi_c <= hadron_fid_phi_max(theta, sector - 1)) return true;

  return false;
}

bool Cuts::Pip(int part) {
  bool _pip = true;
  _pip &= (_data->q(part) == POSITIVE);
  _pip &= Hardon_fid_arjun(part);
  _pip &= dt_Pip_cut(_dt->Get_dt_Pi(part), _data->p(part));
  return _pip;
}
bool Cuts::Pim(int part) {
  bool _pim = true;
  _pim &= (_data->q(part) == NEGATIVE);
  _pim &= Hardon_fid_arjun(part);
  _pim &= dt_Pip_cut(_dt->Get_dt_Pi(part), _data->p(part));
  return _pim;
}
bool Cuts::Prot(int part) {
  bool _prot = true;
  _prot &= (_data->q(part) == POSITIVE);
  _prot &= Hardon_fid_arjun(part);
  _prot &= dt_P_cut(_dt->Get_dt_P(part), _data->p(part));
  return _prot;
}

bool Cuts::Fid_cut() {
  Set_elec_fid();
  return elec_fid_cut();
}

bool Cuts::Beam_cut() {
  bool _beam = true;

  _beam &= (_data->dc_vx(0) > 0.2f && _data->dc_vx(0) < 0.4f);
  _beam &= (_data->dc_vy(0) > -0.1f && _data->dc_vy(0) < 0.16f);
  _beam &= (_data->dc_vz(0) > -2.0f && _data->dc_vz(0) < 2.0f);

  for (short i = 0; i < _data->gpart(); i++) {
    _beam &= (_data->vx(i) > -1.0f && _data->vx(i) < 1.0f);
    _beam &= (_data->vy(i) > -1.5f && _data->vy(i) < 1.5f);
    _beam &= (_data->vz(i) > -2.0f && _data->vz(i) < 2.0f);
  }

  return _beam;
}

bool Cuts::isStrictElecctron() {
  bool _elec = true;
  _elec &= isElecctron();
  // remove CC hit to both
  ////_elec &= (_data->cc_segm(0) / 1000 - 1 != 0);
  // Cut low number of photo electrons in cc
  _elec &= (_data->nphe(0) > 15);

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

bool Cuts::Electron_fid_arjun() {
  float c1e = 12.0;
  float c2e = 18.5;
  float c3e = 0.25;
  float c4e = 15.0;
  float factor_e = 0.416667;
  float p_shift_e = 0.14;

  // Calculates the center phi
  Set_elec_fid();

  // Creating the cut function
  float theta_cut = c1e + c2e / (_data->p(0) + p_shift_e);
  float expon = c3e * pow(_data->p(0), factor_e);
  float del_phi = c4e * pow(sinf((_theta - theta_cut) * DEG2RAD), expon);

  if (abs(_phi_cent) <= del_phi && _theta >= theta_cut) return true;

  return false;
}