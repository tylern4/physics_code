/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef REACTION_H_GUARD
#define REACTION_H_GUARD

#include <iomanip>
#include <iostream>
#include "branches.hpp"
#include "constants.hpp"
#include "physics.hpp"

class Reaction {
 private:
  double _beam_energy = E1D_E0;
  std::unique_ptr<LorentzVector> _beam;
  std::unique_ptr<LorentzVector> _elec;
  std::unique_ptr<LorentzVector> _gamma;
  std::unique_ptr<LorentzVector> _target = std::make_unique<LorentzVector>(0.0, 0.0, 0.0, MASS_P);
  std::unique_ptr<LorentzVector> _prot;
  std::unique_ptr<LorentzVector> _pip;
  std::unique_ptr<LorentzVector> _pim;
  std::unique_ptr<LorentzVector> _neutron;
  std::unique_ptr<LorentzVector> _com;
  std::unique_ptr<LorentzVector> _elec_boosted;
  std::unique_ptr<LorentzVector> _gamma_boosted;
  std::unique_ptr<LorentzVector> _beam_boosted;
  std::unique_ptr<LorentzVector> _pip_boosted;
  std::unique_ptr<LorentzVector> _p_boosted;

  std::shared_ptr<Branches> _data;

  bool _hasE = false;
  bool _hasP = false;
  bool _hasPip = false;
  bool _hasPim = false;
  bool _hasOther = false;
  bool _hasNeutron = false;

  bool _boosted = false;

  short _numProt = 0;
  short _numPip = 0;
  short _numPim = 0;
  short _numPos = 0;
  short _numNeg = 0;
  short _numNeutral = 0;
  short _numOther = 0;

  float _MM = NAN;
  float _MM2 = NAN;

  float _W = NAN;
  float _Q2 = NAN;

  float _theta_e = NAN;
  float _theta_star = NAN;
  float _phi_star = NAN;

  std::string _type = "NAN";

 public:
  Reaction(std::shared_ptr<Branches> data);
  Reaction(std::shared_ptr<Branches> data, bool MC);
  ~Reaction();

  void SetProton(int i);
  void SetPip(int i);
  void SetPim(int i);
  void SetOther(int i);
  void SetNeutron(int i);

  void CalcMissMass();
  double MM();
  double MM2();
  void boost();
  int Type();

  inline double W() { return _W; }
  inline double Q2() { return _Q2; }

  inline double Theta_star() {
    if (!_boosted) boost();
    return _theta_star;
  }
  inline double Phi_star() {
    if (!_boosted) boost();
    return _phi_star;
  }
  inline double Theta_E() {
    if (!_boosted) boost();
    return _theta_e;
  }

  inline bool TwoPion() {
    return ((_numPip == 1 && _numPim == 1) && (_hasE && !_hasP && _hasPip && _hasPim && !_hasNeutron && !_hasOther));
  }
  inline bool ProtonPim() {
    return ((_numProt == 1 && _numPim == 1) && (_hasE && _hasP && !_hasPip && _hasPim && !_hasNeutron && !_hasOther));
  }
  inline bool SinglePip() {
    return ((_numPip == 1) && (_hasE && !_hasP && _hasPip && !_hasPim && !_hasNeutron && !_hasOther));
  }
  inline bool SingleP() {
    return ((_numProt == 1) && (_hasE && _hasP && !_hasPip && !_hasPim && !_hasNeutron && !_hasOther));
  }
  inline bool NeutronPip() {
    return ((_numPip == 1 && _numNeutral == 1) &&
            (_hasE && !_hasP && _hasPip && !_hasPim && _hasNeutron && !_hasOther));
  }

  inline bool channel() { return (this->SinglePip() || this->NeutronPip()) && this->MM_cut(); }

  inline bool MM_cut() {
    bool mm_cut = true;
    // mm_cut &= (this->MM() < 0.987669);
    // mm_cut &= (this->MM() > 0.923374);
    mm_cut &= (this->MM() >= 0.8);
    mm_cut &= (this->MM() <= 1.1);
    return mm_cut;
  }

  inline LorentzVector e_mu() { return *_beam; }
  inline LorentzVector e_mu_prime() { return *_elec; }
  inline LorentzVector p_prime() { return *_prot; }
  inline LorentzVector gamma() { return *_gamma; }
};

#endif
