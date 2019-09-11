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

Reaction::Reaction(std::shared_ptr<Branches> data) : _data(std::move(data)) {
  if (getenv("BEAM_E") != NULL) _beam_energy = atof(getenv("BEAM_E"));
  _beam = std::make_unique<LorentzVector>(0.0, 0.0, _beam_energy, MASS_E);
  _elec = std::make_unique<LorentzVector>(_data->px(0), _data->py(0), _data->pz(0), MASS_E);
  _hasE = true;
  _gamma = std::make_unique<LorentzVector>(0, 0, 0, 0);
  *_gamma = *_beam - *_elec;
  _W = physics::W_calc(*_gamma);
  _Q2 = physics::Q2_calc(*_gamma);
}

Reaction::Reaction(std::shared_ptr<Branches> data, bool MC) : _data(std::move(data)) {
  if (getenv("BEAM_E") != NULL) _beam_energy = atof(getenv("BEAM_E"));
  _beam = std::make_unique<LorentzVector>(0.0, 0.0, _beam_energy, MASS_E);
  _gamma = std::make_unique<LorentzVector>();
  if (MC)
    _elec = std::make_unique<LorentzVector>(_data->pxpart(0), _data->pypart(0), _data->pzpart(0), MASS_E);
  else
    _elec = std::make_unique<LorentzVector>(_data->px(0), _data->py(0), _data->pz(0), MASS_E);

  _hasE = true;
  _gamma = std::make_unique<LorentzVector>(*_beam - *_elec);
  _W = physics::W_calc(*_gamma);
  _Q2 = physics::Q2_calc(*_gamma);
}

Reaction::~Reaction() = default;

void Reaction::SetProton(int i) {
  _numProt++;
  _numPos++;
  _hasP = true;
  _prot = std::make_unique<LorentzVector>(_data->px(i), _data->py(i), _data->pz(i), MASS_P);
}
void Reaction::SetPip(int i) {
  _numPip++;
  _numPos++;
  _hasPip = true;
  _pip = std::make_unique<LorentzVector>(_data->px(i), _data->py(i), _data->pz(i), MASS_PIP);
}
void Reaction::SetPim(int i) {
  _numPim++;
  _numNeg++;
  _hasPim = true;
  _pim = std::make_unique<LorentzVector>(_data->px(i), _data->py(i), _data->pz(i), MASS_PIM);
}

void Reaction::SetNeutron(int i) {
  _numNeutral++;
  _hasNeutron = true;
  _neutron = std::make_unique<LorentzVector>(_data->px(i), _data->py(i), _data->pz(i), MASS_N);
}

void Reaction::SetOther(int i) {
  if (_data->id(i) == NEUTRON)
    Reaction::SetNeutron(i);
  else if (_data->id(i) == PHOTON) {
    _photons.push_back(std::make_unique<LorentzVector>(_data->px(i), _data->py(i), _data->pz(i), 0));
    _numPhotons++;
  } else {
    _numOther++;
    _hasOther = true;
  }
}

void Reaction::CalcMissMass() {
  auto mm = std::make_unique<LorentzVector>();
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
  auto pi0 = std::make_unique<LorentzVector>();
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
  } else if (this->TwoPion()) {
    return 22;
  } else if (this->SingleP()) {
    return 33;
  } else if (this->ProtonPim()) {
    return 44;
  }

  return -1;
}

void Reaction::boost() {
  // Angles gotten from picture in KPark thesis page 11
  // May be wrong still, need to check
  if (_boosted) return;
  _boosted = true;
  if (!(this->SinglePip() || this->NeutronPip())) {
    _theta_e = NAN;
    _theta_star = NAN;
    _phi_star = NAN;
    return;
  }

  auto _com_ = *_target + (*_beam - *_elec);
  auto com = std::make_unique<TLorentzVector>(_com_.X(), _com_.Y(), _com_.Z(), _com_.E());
  auto elec_boosted = std::make_unique<TLorentzVector>(_elec->X(), _elec->Y(), _elec->Z(), _elec->E());
  auto gamma_boosted = std::make_unique<TLorentzVector>(_gamma->X(), _gamma->Y(), _gamma->Z(), _gamma->E());
  auto beam_boosted = std::make_unique<TLorentzVector>(_beam->X(), _beam->Y(), _beam->Z(), _beam->E());
  auto pip_boosted = std::make_unique<TLorentzVector>(_pip->X(), _pip->Y(), _pip->Z(), _pip->E());

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

  auto _temp = std::make_unique<LorentzVector>(pip_boosted->X(), pip_boosted->Y(), pip_boosted->Z(), pip_boosted->M());
  _theta_e = elec_boosted->Theta();
  _theta_star = _temp->Theta();
  _phi_star = physics::phi_boosted(_temp);
}
