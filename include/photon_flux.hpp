/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#ifndef PF_H
#define PF_H
#include <iostream>
#include "color.hpp"
#include "constants.hpp"
#include "physics.hpp"

class PhotonFlux {
 private:
  double _beam_energy = E1D_E0;
  double _target_mass = MASS_P;
  double _beam_momentum;
  double _nu;
  double _scattered_energy;
  double _scattered_momentum;
  double _W;
  double _Q2;
  double _flux;

  double Momentum(double E, double M);
  double photon_energy();
  double theta_calc();
  double epsilon_calc();
  double photon_flux();

 public:
  PhotonFlux();
  PhotonFlux(double W, double Q2);
  PhotonFlux(LorentzVector e_mu, LorentzVector e_mu_prime);
  ~PhotonFlux();

  double GetVirtualPhotonFlux();
  std::string FluxString();
  void PrintFlux();
};

#endif
