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

  //_phi_cent = _phi + phi_center[_sec];

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

bool Cuts::Sanity() {
  bool _elec = true;
  _elec &= (_data->gpart() > 0);  // Number of good particles is greater than 0
  _elec &= (_data->gpart() < 5);
  if (!_elec) return false;
  _elec &= (_data->q(0) == NEGATIVE);
  _elec &= (_data->ec(0) > 0);
  _elec &= (_data->cc(0) > 0);
  _elec &= (_data->stat(0) > 0);  // First Particle stat
  _elec &= (_data->sc(0) > 0);
  _elec &= (_data->dc(0) > 0);
  _elec &= (_data->dc_stat(0) > 0);

  return _elec;
}

bool Cuts::isElecctron() {
  bool _elec = true;
  _elec &= (_data->gpart() > 0);  // Number of good particles is greater than 0
  _elec &= (_data->gpart() < 5);
  if (!_elec) return false;
  _elec &= (_data->q(0) == NEGATIVE);

  _elec &= (_data->ec(0) > 0);
  _elec &= (_data->cc(0) > 0);
  _elec &= (_data->stat(0) > 0);  // First Particle stat
  _elec &= (_data->sc(0) > 0);
  _elec &= (_data->dc(0) > 0);
  _elec &= (_data->dc_stat(0) > 0);
  //_elec &= _data->nphe(0) < 20;

  // Cut out low ec inner
  _elec &= (_data->ec_ei(0) >= 0.05);
  // Minimum momentum cut
  _elec &= (_data->p(0) > MIN_P_CUT);

  return _elec;
}

bool e1d_Cuts::isElecctron() {
  bool _elec = true;
  _elec &= Cuts::isElecctron();
  _elec &= e1d_Cuts::Beam_cut();
  // Sampling fraction cut
  _elec &= e1d_Cuts::sf_cut(_data->etot(0) / _data->p(0), _data->p(0));
  if (!_elec) return _elec;
  // Fid Cuts
  _elec &= fid_chern_cut();
  if (!_elec) return _elec;
  _elec &= Fid_cut();

  return _elec;
}

bool e1f_Cuts::isElecctron() {
  bool _elec = true;
  _elec &= Cuts::isElecctron();
  _elec &= e1f_Cuts::Beam_cut();
  //_elec &= e1f_Cuts::sf_cut(_data->etot(0) / _data->p(0), _data->p(0));
  if (!_elec) return _elec;
  // Fid Cuts
  _elec &= fid_chern_cut();
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
  _pip &= dt_Pip_cut(part);
  return _pip;
}
bool Cuts::Pim(int part) {
  bool _pim = true;
  _pim &= (_data->q(part) == NEGATIVE);
  _pim &= Hardon_fid_arjun(part);
  _pim &= dt_Pip_cut(part);
  return _pim;
}
bool Cuts::Prot(int part) {
  bool _prot = true;
  _prot &= (_data->q(part) == POSITIVE);
  _prot &= Hardon_fid_arjun(part);
  _prot &= dt_P_cut(part);
  return _prot;
}

bool Cuts::Pipish(int part) {
  bool _pip = true;
  _pip &= (_data->q(part) == POSITIVE);
  _pip &= Hardon_fid_arjun(part);
  return _pip;
}
bool Cuts::Pimish(int part) {
  bool _pim = true;
  _pim &= (_data->q(part) == NEGATIVE);
  _pim &= Hardon_fid_arjun(part);
  return _pim;
}
bool Cuts::Protish(int part) {
  bool _prot = true;
  _prot &= (_data->q(part) == POSITIVE);
  _prot &= Hardon_fid_arjun(part);
  return _prot;
}

bool Cuts::Fid_cut() {
  Set_elec_fid();
  return elec_fid_cut();
}

bool Cuts::Beam_cut() { return true; }

bool e1d_Cuts::Beam_cut() {
  bool _beam = true;

  _beam &= (_data->dc_vx(0) > 0.2f && _data->dc_vx(0) < 0.4f);
  _beam &= (_data->dc_vy(0) > -0.1f && _data->dc_vy(0) < 0.16f);
  _beam &= (_data->dc_vz(0) > -5.0f && _data->dc_vz(0) < 5.0f);

  for (short i = 0; i < _data->gpart(); i++) {
    _beam &= (_data->vx(i) > -1.0f && _data->vx(i) < 1.0f);
    _beam &= (_data->vy(i) > -1.5f && _data->vy(i) < 1.5f);
    _beam &= (_data->vz(i) > -5.0f && _data->vz(i) < 5.0f);
  }

  return _beam;
}

bool e1f_Cuts::Beam_cut() {
  bool _beam = true;

  _beam &= (abs(_data->dc_vx(0)) < 0.3f);
  _beam &= (abs(_data->dc_vy(0)) < 0.4f);
  //_beam &= (_data->dc_vz(0) > -3.0f && _data->dc_vz(0) < 0.0f);
  /*
    for (short i = 0; i < _data->gpart(); i++) {
      _beam &= (abs(_data->vx(i)) < 0.2f);
      _beam &= (abs(_data->vy(i)) < 0.4f);
      _beam &= (_data->vz(i) > -3.0f && _data->vz(i) < 0.0f);
    }
  */
  return _beam;
}

bool Cuts::isStrictElecctron() {
  bool _elec = true;
  _elec &= isElecctron();
  _elec &= (_data->nphe(0) > 15);
  return _elec;
}

float e1d_Cuts::sf_top_fit(float P) {
  double par[3] = {0.368209, 0.000961273, 4.8e-07};
  // double par[3] = {0.3269, 0.000336, 7.731e-7};
  double x[1] = {P};
  return func::ec_fit_func(x, par);
}
float e1d_Cuts::sf_bot_fit(float P) {
  double par[3] = {0.162189, 0.0134756, -2e-05};
  // double par[3] = {0.1787, 0.02032, -2.726e-6};
  double x[1] = {P};
  return func::ec_fit_func(x, par);
}
bool e1d_Cuts::sf_cut(float sf, float P) { return ((sf > e1d_Cuts::sf_bot_fit(P)) && (sf < e1d_Cuts::sf_top_fit(P))); }

float e1f_Cuts::sf_top_fit(float P) {
  double par[3] = {0.405, 0.007693, 4.993e-07};
  double x[1] = {P};
  return func::ec_fit_func(x, par);

  return true;
}
float e1f_Cuts::sf_bot_fit(float P) {
  double par[3] = {0.1272, 0.04403, -1.998e-05};
  double x[1] = {P};
  return func::ec_fit_func(x, par);

  return true;
}
bool e1f_Cuts::sf_cut(float sf, float P) { return ((sf > e1f_Cuts::sf_bot_fit(P)) && (sf < e1f_Cuts::sf_top_fit(P))); }

bool Cuts::dt_P_cut(int i) {
  float dt = _dt->Get_dt_P(i);
  int sec = _data->dc_sect(i) - 1;
  if (sec == -1) return false;
  float p = _data->p(i);
  bool _cut = true;
  //_cut &= (dt <= func::dt_poly4(dt_P_const_top, p));
  //_cut &= (dt >= func::dt_poly4(dt_P_const_bottom, p));

  _cut &= (dt <= func::log_sqrt_pol1(dt_pip_top, p));
  _cut &= (dt >= func::log_sqrt_pol1(dt_pip_bottom, p));

  return _cut;
}

bool Cuts::dt_K_cut(int i) {
  float dt = _dt->Get_dt_K(i);
  int sec = _data->dc_sect(i) - 1;
  if (sec == -1) return false;
  float p = _data->p(i);
  bool _cut = true;
  _cut &= (dt <= func::dt_poly4(dt_pip_const_top, p));
  _cut &= (dt >= func::dt_poly4(dt_pip_const_bottom, p));

  return _cut;
}

bool Cuts::dt_Pip_cut(int i) {
  float dt = _dt->Get_dt_Pi(i);
  short sec = _data->dc_sect(i) - 1;
  if (sec == -1) return false;
  float p = _data->p(i);
  bool _cut = true;

  //_cut &= (dt <= func::dt_poly4(dt_pip_const_top, p));
  //_cut &= (dt >= func::dt_poly4(dt_pip_const_bottom, p));

  _cut &= (dt <= func::log_sqrt_pol1(dt_pip_top, p));
  _cut &= (dt >= func::log_sqrt_pol1(dt_pip_bottom, p));

  return _cut;
}

bool Cuts::elec_fid_cut() {
  // float c[3] = {0.0548311203, 0.0327878012, 17.0683287374};
  float c[3] = {0.04, 0.03, 15};
  float y = c[0] * _phi_cent * _phi_cent + c[1] * _phi_cent + c[2];

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

bool Cuts::fid_chern_cut() {
  float A = -0.000785;
  float B = 0;
  float C = -0.00168;
  float D = 1;

  auto p0_vec = TVector3(_data->dc_xsc(0), _data->dc_ysc(0), _data->dc_zsc(0));
  auto n_vec = TVector3(_data->dc_cxsc(0), _data->dc_cysc(0), _data->dc_czsc(0));
  auto S_vec = TVector3(A, B, C);

  auto numer = A * _data->dc_xsc(0) + B * _data->dc_ysc(0) + C * _data->dc_zsc(0) + D;
  auto denom = S_vec.Dot(n_vec);

  auto t_vec = n_vec * abs(numer / denom);

  p0_vec += t_vec;
  float _cc_theta = acosf(p0_vec.Z() / p0_vec.Mag());
  float _cc_phi = atanf(p0_vec.Y() / p0_vec.X());

  float _cc_x = _data->cc_r(0) * sinf(_cc_theta) * cosf(_cc_phi);
  float _cc_y = _data->cc_r(0) * sinf(_cc_theta) * sinf(_cc_phi);
  return func::fid_chern(_cc_x, _cc_y);
}