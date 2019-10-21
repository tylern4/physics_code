/**************************************/
/*                                    */
/*  Created by Nick Tyler             */
/*	University Of South Carolina  */
/**************************************/
#include "reaction.hpp"
#include <mutex>
#include "TLorentzRotation.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include "TRandom.h"

Reaction::Reaction(const std::shared_ptr<Branches>& data) : _data(data) {
  _hasE = true;
  _sector = data->dc_sect(0);
  _elec = physics::fourVec(_data->px(0), _data->py(0), _data->pz(0), MASS_E);

  // this->correct_mom();
  *_gamma = *_beam - *_elec;
  _W = physics::W_calc(*_gamma);
  _Q2 = physics::Q2_calc(*_gamma);
  _xb = physics::xb_calc(*_gamma);
}

Reaction::~Reaction() = default;

void Reaction::correct_mom() {
  for (int i = 0; i < 6; i++)
    for (int ii = 0; ii < 16; ii++) par[i][ii] = 0.0;
  float torus = 1.0;
  float phi = _elec->Phi() * RAD2DEG;
  // transforming phi: [-30,270] -> [-30,30]
  if (phi < -30.) phi = phi + 360.;
  phi = ((phi + 30.) / 60. - (int(phi + 30.)) / 60) * 60. - 30.;
  phi = phi * DEG2RAD;

  float factor = _elec->Theta() / (sin(4. * _elec->Theta()) * sin(4. * _elec->Theta()));
  if (_elec->Theta() * RAD2DEG >= 22.5) factor = _elec->Theta();

  float theta = _elec->Theta();

  float dtheta = (par[_sector - 1][0] + par[_sector - 1][1] * phi) * cos(theta) / cos(phi) +
                 (par[_sector - 1][2] + par[_sector - 1][3] * phi) * sin(theta) +
                 par[_sector - 1][14] * cos(theta) * cos(theta) / cos(phi);

  // new corrected theta
  theta += dtheta;

  float term = (par[_sector - 1][4] + par[_sector - 1][5] * phi) * cos(theta) / cos(phi) +
               (par[_sector - 1][6] + par[_sector - 1][7] * phi) * sin(theta) +
               par[_sector - 1][15] * cos(theta) * cos(theta) / cos(phi);

  term *= _elec->P() * 3375. / (_data->q(0) * torus * 2250.) * factor;
  float dmom =
      term + par[_sector - 1][8] * cos(theta) + par[_sector - 1][9] * sin(theta) +
      par[_sector - 1][10] * sin(2. * theta) +
      (par[_sector - 1][11] * cos(theta) + par[_sector - 1][12] * sin(theta) + par[_sector - 1][13] * sin(2. * theta)) *
          phi;
  // new corrected momentum
  float p_cor = _elec->P() * (1. + dmom);
  phi = _elec->Phi();

  //_elec->SetPxPyPzE(p_cor * cos(phi) * sin(theta), p_cor * sin(phi) * sin(theta), p_cor * cos(theta), sqrt(p_cor *
  // p_cor + MASS_E * MASS_E));
}

void Reaction::SetProton(int i) {
  _numProt++;
  _numPos++;
  _hasP = true;
  _prot = physics::fourVec(_data->px(i), _data->py(i), _data->pz(i), MASS_P);
}
void Reaction::SetPip(int i) {
  _numPip++;
  _numPos++;
  _hasPip = true;
  _pip = physics::fourVec(_data->px(i), _data->py(i), _data->pz(i), MASS_PIP);
}
void Reaction::SetPim(int i) {
  _numPim++;
  _numNeg++;
  _hasPim = true;
  _pim = physics::fourVec(_data->px(i), _data->py(i), _data->pz(i), MASS_PIM);
}

void Reaction::SetNeutron(int i) {
  _numNeutral++;
  _hasNeutron = true;
  _neutron = physics::fourVec(_data->px(i), _data->py(i), _data->pz(i), MASS_N);
}

void Reaction::SetOther(int i) {
  if (_data->id(i) == NEUTRON)
    Reaction::SetNeutron(i);
  else if (_data->id(i) == PHOTON) {
    _photons.push_back(physics::fourVec(_data->px(i), _data->py(i), _data->pz(i), 0));
    _numPhotons++;
  } else {
    _numOther++;
    _hasOther = true;
  }
}

void Reaction::CalcMissMass() {
  auto mm = std::make_shared<LorentzVector>();
  *mm += (*_beam - *_elec + *_target);
  if (SinglePip() || NeutronPip()) {
    *mm -= *_pip;
    _MM = mm->M();
    _MM2 = mm->M2();
  } else if (TwoPion()) {
    *mm -= *_pip;
    *mm -= *_pim;
    _MM = mm->mag();
    _MM2 = mm->mag2();
  } else if (ProtonPim()) {
    *mm -= *_prot;
    *mm -= *_pim;
    _MM = mm->mag();
    _MM2 = mm->mag2();
  } else if (SingleP()) {
    *mm -= *_prot;
    _MM = mm->mag();
    _MM2 = mm->mag2();
  }
  auto pi0 = std::make_shared<LorentzVector>();
  if (_numPhotons >= 2) {
    for (auto& p : _photons) {
      *pi0 += *p;
    }
    _pi0_mass = pi0->mag();
    _pi0_mass2 = pi0->mag2();
  }
}

double Reaction::MM() {
  if (!_MM_calc) {
    CalcMissMass();
    _MM_calc = true;
  }
  return _MM;
}
double Reaction::MM2() {
  if (!_MM_calc) {
    CalcMissMass();
    _MM_calc = true;
  }
  return _MM2;
}

float Reaction::pi0_mass() {
  if (!_MM_calc) {
    CalcMissMass();
    _MM_calc = true;
  }
  return _pi0_mass;
}

float Reaction::pi0_mass2() {
  if (!_MM_calc) {
    CalcMissMass();
    _MM_calc = true;
  }
  return _pi0_mass2;
}

int Reaction::Type() {
  if (this->SinglePip()) {
    return 0;
  } else if (this->NeutronPip()) {
    return 10;
  } else if (this->SingleP()) {
    return 22;
  } else if (this->PPi0()) {
    return 222;
  } else if (this->ProtonPim()) {
    return 3333;
  }

  return -1;
}

void Reaction::boost() {
  // Angles gotten from picture in KPark thesis page 11
  // May be wrong still, need to check
  if (_boosted) return;
  _boosted = true;

  if (!(this->SinglePip() || this->NeutronPip() || this->SingleP() || this->PPi0())) {
    _theta_e = NAN;
    _theta_star = NAN;
    _phi_star = NAN;
    return;
  }

  auto _com_ = *_target + (*_beam - *_elec);
  auto com = std::make_shared<TLorentzVector>(_com_.X(), _com_.Y(), _com_.Z(), _com_.E());
  auto elec_boosted = std::make_shared<TLorentzVector>(_elec->X(), _elec->Y(), _elec->Z(), _elec->E());
  auto gamma_boosted = std::make_shared<TLorentzVector>(_gamma->X(), _gamma->Y(), _gamma->Z(), _gamma->E());
  auto beam_boosted = std::make_shared<TLorentzVector>(_beam->X(), _beam->Y(), _beam->Z(), _beam->E());
  auto pip_boosted = std::make_shared<TLorentzVector>(_pip->X(), _pip->Y(), _pip->Z(), _pip->E());

  //! Varsets
  //! Calculate rotation: taken from Evan's phys-ana-omega on 08-05-13
  // Copied and modified from Arjun's code
  TVector3 uz = gamma_boosted->Vect().Unit();
  TVector3 ux = (beam_boosted->Vect().Cross(elec_boosted->Vect())).Unit();
  ux.Rotate(-PI / 2, uz);
  TRotation r3;  // = new TRotation();
  r3.SetZAxis(uz, ux).Invert();
  //! _w and _q are in z-direction
  TVector3 boost(-1 * com->BoostVector());
  TLorentzRotation r4(r3);  //*_boost);
  r4 *= boost;              //*_3rot;

  gamma_boosted->Transform(r4);
  elec_boosted->Transform(r4);
  beam_boosted->Transform(r4);
  pip_boosted->Transform(r4);

  auto _temp = physics::fourVec(pip_boosted->X(), pip_boosted->Y(), pip_boosted->Z(), pip_boosted->M());
  _theta_e = elec_boosted->Theta();
  _theta_star = _temp->Theta();
  _phi_star = physics::phi_boosted(_temp);
}

void Reaction::_boost() {
  // Angles gotten from picture in KPark thesis page 11
  // May be wrong still, need to check
  if (_boosted) return;
  _boosted = true;

  if (!(this->SinglePip() || this->NeutronPip() || this->SingleP() || this->PPi0())) {
    _theta_e = NAN;
    _theta_star = NAN;
    _phi_star = NAN;
    return;
  }

  auto com = *_target + (*_beam - *_elec);
  auto uz = _gamma->BoostToCM();
  auto ux = (_beam->Vect().Cross(_elec->Vect())).Unit();
  auto boost = (-1 * com.BoostToCM());

  //_elec

  _theta_e = _elec->Theta();
  _theta_star = _pip->Theta();
  _phi_star = physics::phi_boosted(_pip);
}

MCReaction::MCReaction(std::shared_ptr<Branches> data) : Reaction(data) {
  _elec_thrown = physics::fourVec(_data->pxpart(0), _data->pypart(0), _data->pzpart(0), MASS_E);
  _gamma_thrown = std::make_shared<LorentzVector>(*_beam - *_elec_thrown);
  _W_thrown = physics::W_calc(*_gamma_thrown);
  _Q2_thrown = physics::Q2_calc(*_gamma_thrown);
  _pip_thrown = physics::fourVec(_data->pxpart(1), _data->pypart(1), _data->pzpart(1), MASS_PIP);
}

float MCReaction::Theta_E() { return _elec_thrown->Theta(); }
float MCReaction::Theta_star() { return _pip_thrown->Theta(); }
float MCReaction::Phi_star() { return physics::phi_boosted(_pip_thrown); }
