/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "reaction.hpp"

Reaction::Reaction() {
  _beam = std::make_unique<TLorentzVector>();
  if (getenv("BEAM_E") != NULL) _beam_energy = atof(getenv("BEAM_E"));
  _beam->SetPxPyPzE(0.0, 0.0, sqrt(_beam_energy * _beam_energy - MASS_E * MASS_E), _beam_energy);

  _gamma = std::make_unique<TLorentzVector>();
  _target = std::make_unique<TLorentzVector>(0.0, 0.0, 0.0, MASS_P);
  _elec = std::make_unique<TLorentzVector>();
  _prot = std::make_unique<TLorentzVector>();
  _pip = std::make_unique<TLorentzVector>();
  _pim = std::make_unique<TLorentzVector>();
}

Reaction::~Reaction() {}

void Reaction::SetElec(float p, float cx, float cy, float cz) {
  _hasE = true;
  _elec->SetXYZM(p * cx, p * cy, p * cz, MASS_E);

  *_gamma += *_beam - *_elec;

  // Can calculate W and Q2 here
  _W = physics::W_calc(*_beam, *_elec);
  _Q2 = physics::Q2_calc(*_beam, *_elec);
}

void Reaction::SetProton(float p, float cx, float cy, float cz) {
  _numP++;
  _numPos++;
  _hasP = true;
  _prot->SetXYZM(p * cx, p * cy, p * cz, MASS_P);
}
void Reaction::SetPip(float p, float cx, float cy, float cz) {
  _numPip++;
  _numPos++;
  _hasPip = true;
  _pip->SetXYZM(p * cx, p * cy, p * cz, MASS_PIP);
}
void Reaction::SetPim(float p, float cx, float cy, float cz) {
  _numPim++;
  _numNeg++;
  _hasPim = true;
  _pim->SetXYZM(p * cx, p * cy, p * cz, MASS_PIM);
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

float Reaction::MM() {
  if (_MM != _MM) CalcMissMass();
  return _MM;
}
float Reaction::MM2() {
  if (_MM2 != _MM2) CalcMissMass();
  return _MM2;
}
