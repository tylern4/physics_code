/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef REACTION_H_GUARD
#define REACTION_H_GUARD

#include <iostream>
#include <map>
#include "TLorentzVector.h"
#include "branches.hpp"
#include "constants.hpp"
#include "particle.hpp"
#include "physics.hpp"

class Reaction {
 private:
  std::map<int, double> _mass_map = {{PROTON, MASS_P}, {-PROTON, MASS_P},  {NEUTRON, MASS_N},  {PIP, MASS_PIP},
                                     {PIM, MASS_PIM},  {PI0, MASS_PI0},    {KP, MASS_KP},      {KM, MASS_KM},
                                     {PHOTON, MASS_G}, {ELECTRON, MASS_E}, {-ELECTRON, MASS_E}};
  double _beam_energy = E1D_E0;
  std::unique_ptr<TLorentzVector> _beam;
  std::unique_ptr<TLorentzVector> _elec;
  std::unique_ptr<TLorentzVector> _gamma;
  std::unique_ptr<TLorentzVector> _target;
  std::unique_ptr<TLorentzVector> _prot;
  std::unique_ptr<TLorentzVector> _pip;
  std::unique_ptr<TLorentzVector> _pim;
  std::unique_ptr<TLorentzVector> _other;
  std::unique_ptr<TLorentzVector> _neutron;

  /*
    std::unique_ptr<Particle> _beam_particle;
    std::unique_ptr<Particle> _elec_particle;
    std::unique_ptr<Particle> _gamma_particle;
    std::unique_ptr<Particle> _target_particle;
    std::unique_ptr<Particle> _prot_particle;
    std::unique_ptr<Particle> _pip_particle;
    std::unique_ptr<Particle> _pim_particle;
    std::unique_ptr<Particle> _other_particle;
    std::unique_ptr<Particle> _neutron_particle;
  */
  Branches *_data;

  bool _hasE = false;
  bool _hasP = false;
  bool _hasPip = false;
  bool _hasPim = false;
  bool _hasOther = false;
  bool _hasNeutron = false;

  short _numProt = 0;
  short _numPip = 0;
  short _numPim = 0;
  short _numPos = 0;
  short _numNeg = 0;
  short _numNeutral = 0;
  short _numOther = 0;

  float _MM = std::nanf("-99");
  float _MM2 = std::nanf("-99");

  float _W = std::nanf("-99");
  float _Q2 = std::nanf("-99");

 public:
  Reaction(Branches *data);
  ~Reaction();

  void SetProton(int i);
  void SetPip(int i);
  void SetPim(int i);
  void SetOther(int i);
  void SetNeutron(int i);

  void CalcMissMass();
  float MM();
  float MM2();

  inline float W() { return _W; }
  inline float Q2() { return _Q2; }

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

  inline TLorentzVector e_mu() { return *_beam; }
  inline TLorentzVector e_mu_prime() { return *_elec; }
  inline TLorentzVector gamma() { return *_gamma; }
};

#endif
