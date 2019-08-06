/**************************************/
/*                                    */
/*  Created by Nick Tyler             */
/*	University Of South Carolina  */
/**************************************/
#include "reaction.hpp"

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
  _elec = std::make_unique<LorentzVector>(_data->pxpart(0), _data->pypart(0), _data->pzpart(0), MASS_E);
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
    SetNeutron(i);
  else {
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
}

double Reaction::MM() {
  // if (std::isnan(_MM))
  CalcMissMass();
  return _MM;
}
double Reaction::MM2() {
  // if (std::isnan(_MM2))
  CalcMissMass();
  return _MM2;
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
  _boosted = true;
  /*
  _com = std::make_unique<LorentzVector>(*_target);
  *_com += (*_beam - *_elec);
  _elec_boosted = std::make_unique<LorentzVector>(*_elec);
  _gamma_boosted = std::make_unique<LorentzVector>(*_gamma);
  _beam_boosted = std::make_unique<LorentzVector>(*_beam);
  _pip_boosted = std::make_unique<LorentzVector>(*_pip);
  _p_boosted = std::make_unique<LorentzVector>(*_prot);

  //! Varsets
  //! Calculate rotation: taken from Evan's phys-ana-omega on 08-05-13
  // Copied and modified from Arjun's code
  TVector3 uz = _gamma->Vect().Unit();
  TVector3 ux = (_beam->Vect().Cross(_elec->Vect())).Unit();
  ux.Rotate(-PI / 2, uz);
  TRotation r3;  // = new TRotation();
  r3.SetZAxis(uz, ux).Invert();
  //! _w and _q are in z-direction
  TVector3 boost(-1 * _com->BoostVector());
  TLorentzRotation r4(r3);  //*_boost);
  r4 *= boost;              //*_3rot;

  _gamma_boosted->Transform(r4);
  _elec_boosted->Transform(r4);
  _beam_boosted->Transform(r4);
  _pip_boosted->Transform(r4);
  _p_boosted->Transform(r4);

  if ((this->SinglePip() || this->NeutronPip())) {
    _theta_e = _elec_boosted->Theta() / D2R;
    _theta_star = _pip_boosted->Theta() / D2R;
    _phi_star = physics::phi_boosted(_pip_boosted) / D2R;
  } else if (this->SingleP()) {
    _theta_e = _elec_boosted->Theta() / D2R;
    _theta_star = _p_boosted->Theta() / D2R;
    _phi_star = physics::phi_boosted(_p_boosted) / D2R;
  } else {
    _theta_e = _elec_boosted->Theta();
    _theta_star = NAN;
    _phi_star = NAN;
  }

  */
}
