/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "reaction.hpp"

Reaction::Reaction(std::shared_ptr<Branches> data) : _data(data) {
  _beam = std::make_unique<TLorentzVector>();
  if (getenv("BEAM_E") != NULL) _beam_energy = atof(getenv("BEAM_E"));
  _beam->SetPxPyPzE(0.0, 0.0, sqrt(_beam_energy * _beam_energy - MASS_E * MASS_E), _beam_energy);

  _gamma = std::make_unique<TLorentzVector>();
  _com = std::make_unique<TLorentzVector>();
  _elec_com = std::make_unique<TLorentzVector>();
  _target = std::make_unique<TLorentzVector>(0.0, 0.0, 0.0, MASS_P);
  _elec = std::make_unique<TLorentzVector>();
  _prot = std::make_unique<TLorentzVector>();
  _pip = std::make_unique<TLorentzVector>();
  _pim = std::make_unique<TLorentzVector>();
  //_other = std::make_unique<TLorentzVector>();
  _neutron = std::make_unique<TLorentzVector>();

  _hasE = true;
  _elec->SetXYZM(_data->px(0), _data->py(0), _data->pz(0), MASS_E);
  *_gamma += *_beam - *_elec;
  _W = physics::W_calc(*_gamma);
  _Q2 = physics::Q2_calc(*_gamma);

  *_com += (*_beam - *_elec) + *_target;
  _elec_com->SetXYZM(_data->px(0), _data->py(0), _data->pz(0), MASS_E);
  _elec_com->Boost(-_com->BoostVector());
  _theta_star = _elec_com->Theta() * TMath::RadToDeg();
  _phi_star = _elec_com->Phi() * TMath::RadToDeg();
  //_phi_star = _phi_star < -30 ? _phi_star + 360 : _phi_star;
}

Reaction::~Reaction() {}

void Reaction::SetProton(int i) {
  _numProt++;
  _numPos++;
  _hasP = true;
  _prot->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_P);
}
void Reaction::SetPip(int i) {
  _numPip++;
  _numPos++;
  _hasPip = true;
  _pip->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIP);
}
void Reaction::SetPim(int i) {
  _numPim++;
  _numNeg++;
  _hasPim = true;
  _pim->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIM);
}

void Reaction::SetNeutron(int i) {
  _numNeutral++;
  _hasNeutron = true;
  _neutron->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_N);
}

void Reaction::SetOther(int i) {
  if (_data->id(i) == NEUTRON)
    SetNeutron(i);
  else {
    _numOther++;
    _hasOther = true;
    //_other->push_back(std::make_unique<TLorentzVector>());
    //_other->at(_other->size() - 1)->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), _mass_map[_data->id(i)]);
  }
}

void Reaction::CalcMissMass() {
  auto mm = std::make_unique<TLorentzVector>();
  *mm += (*_beam - *_elec + *_target);
  if (SinglePip()) {
    *mm -= *_pip;
    _MM = mm->M();
    _MM2 = mm->M2();
  } else if (TwoPion()) {
    *mm -= *_prot;
    *mm -= *_pip;
    *mm -= *_pim;
    _MM = mm->M();
    _MM2 = mm->M2();
  } else if (ProtonPim()) {
    *mm -= *_prot;
    *mm -= *_pim;
    _MM = mm->M();
    _MM2 = mm->M2();
  } else if (SingleP()) {
    *mm -= *_prot;
    _MM = mm->M();
    _MM2 = mm->M2();
  }
}

double Reaction::MM() {
  if (_MM != _MM) CalcMissMass();
  return _MM;
}
double Reaction::MM2() {
  if (_MM2 != _MM2) CalcMissMass();
  return _MM2;
}
