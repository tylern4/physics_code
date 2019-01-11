/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef REACTION_H_GUARD
#define REACTION_H_GUARD

#include <iostream>
#include "TLorentzVector.h"
#include "constants.hpp"
#include "physics.hpp"

class Reaction {
 private:
  double _beam_energy = E1D_E0;
  std::unique_ptr<TLorentzVector> _beam;
  std::unique_ptr<TLorentzVector> _elec;
  std::unique_ptr<TLorentzVector> _gamma;
  std::unique_ptr<TLorentzVector> _target;
  std::unique_ptr<TLorentzVector> _prot;
  std::unique_ptr<TLorentzVector> _pip;
  std::unique_ptr<TLorentzVector> _pim;

  bool _hasE = false;
  bool _hasP = false;
  bool _hasPip = false;
  bool _hasPim = false;

  short _numP = 0;
  short _numPip = 0;
  short _numPim = 0;
  short _numPos = 0;
  short _numNeg = 0;

  float _MM = std::nanf("-99");
  float _MM2 = std::nanf("-99");

  float _W = std::nanf("-99");
  float _Q2 = std::nanf("-99");

 public:
  Reaction();
  Reaction(float p, float cx, float cy, float cz);
  ~Reaction();

  void SetElec(float p, float cx, float cy, float cz);
  void SetProton(float p, float cx, float cy, float cz);
  void SetPip(float p, float cx, float cy, float cz);
  void SetPim(float p, float cx, float cy, float cz);

  void CalcMissMass();
  float MM();
  float MM2();

  inline float W() { return _W; }
  inline float Q2() { return _Q2; }

  inline bool TwoPion() { return (_hasE && _hasP && _hasPip && _hasPim); }
  inline bool ProtonPim() { return (_hasE && _hasP && _hasPim && !_hasPip); }
  inline bool SinglePip() { return (_hasE && _hasPip && !_hasP && !_hasPim); }
  inline bool SingleP() { return (_hasE && _hasP && !_hasPip && !_hasPim); }
  inline bool Single_pi() { return ((_numPip + _numPim) == 1 && _numP == 0); }

  inline TLorentzVector e_mu() { return *_beam; }
  inline TLorentzVector e_mu_prime() { return *_elec; }
  inline TLorentzVector gamma() { return *_gamma; }
};

#endif
