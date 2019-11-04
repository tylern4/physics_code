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
 protected:
  double _beam_energy = E1D_E0;
  std::shared_ptr<LorentzVector> _beam = physics::fourVec(0.0, 0.0, BEAM_E, MASS_E);
  std::shared_ptr<LorentzVector> _target = physics::fourVec(0.0, 0.0, 0.0, MASS_P);
  std::shared_ptr<LorentzVector> _gamma = physics::fourVec(0, 0, 0, 0);

  std::shared_ptr<LorentzVector> _elec = nullptr;
  std::shared_ptr<LorentzVector> _prot = nullptr;
  std::shared_ptr<LorentzVector> _pip = nullptr;
  std::shared_ptr<LorentzVector> _pim = nullptr;
  std::shared_ptr<LorentzVector> _neutron = nullptr;
  std::vector<std::shared_ptr<LorentzVector>> _photons;

  std::shared_ptr<Branches> _data;

  double par[6][16];

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
  short _numPhotons = 0;
  short _numOther = 0;

  short _sector = -1;

  bool _MM_calc = false;
  float _MM = NAN;
  float _MM2 = NAN;

  float _pi0_mass = NAN;
  float _pi0_mass2 = NAN;

  float _W = NAN;
  float _Q2 = NAN;
  float _xb = NAN;

  float _theta_e = NAN;
  float _theta_star = NAN;
  float _phi_star = NAN;

  std::string _type = "NAN";

 public:
  Reaction() : _data(nullptr){};
  Reaction(const std::shared_ptr<Branches>& data);

  ~Reaction();
  void correct_mom();
  void SetProton(int i);
  void SetPip(int i);
  void SetPim(int i);
  void SetOther(int i);
  void SetNeutron(int i);

  void CalcMissMass();
  double MM();
  double MM2();
  float pi0_mass();
  float pi0_mass2();
  void boost();
  void _boost_pip();
  void _boost_p();
  int Type();

  friend std::ostream& operator<<(std::ostream& os, Reaction& e) {
    return os << '{' << "W:" << e._W << " Q2:" << e._Q2 << " type:" << e.Type() << '}';
  }

  inline float E_prime() { return _elec->E(); }
  inline short sector() { return _sector; }
  inline float W() { return _W; }
  inline float Q2() { return _Q2; }
  inline float xb() { return _xb; }

  inline float Theta_star() {
    if (!_boosted) boost();
    return _theta_star;
  }
  inline float Phi_star() {
    if (!_boosted) boost();
    return _phi_star;
  }
  inline float Theta_E() {
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

  inline bool PPi0() {
    return (SingleP() &&
            ((this->MM() >= 0.05 && this->MM() <= 0.3) || (this->pi0_mass() >= 0.05 && this->pi0_mass() <= 0.2)));
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

  inline LorentzVector& e_mu() { return *_beam; }
  inline LorentzVector& e_mu_prime() { return *_elec; }
  inline LorentzVector& p_prime() { return *_prot; }
  inline LorentzVector& gamma() { return *_gamma; }
};

class MCReaction : public Reaction {
  float _W_thrown = NAN;
  float _Q2_thrown = NAN;
  std::shared_ptr<LorentzVector> _elec_thrown;
  std::shared_ptr<LorentzVector> _gamma_thrown;
  std::shared_ptr<LorentzVector> _pip_thrown;
  std::shared_ptr<LorentzVector> _neutron_thrown;

 public:
  MCReaction(std::shared_ptr<Branches> data);
  inline float W_thrown() { return _W_thrown; }
  inline float Q2_thrown() { return _Q2_thrown; }

  float Theta_star();
  float Phi_star();
  float Theta_E();
};
#endif
