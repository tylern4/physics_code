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
  float _beam_energy = NAN;
  float _target_mass = MASS_P;
  float _beam_momentum;
  float _nu;
  float _scattered_energy;
  float _scattered_momentum;
  float _W;
  float _Q2;
  float _flux;

  float Momentum(float E, float M);
  float photon_energy();
  float theta_calc();
  float epsilon_calc();
  float photon_flux();

 public:
  PhotonFlux();
  PhotonFlux(float W, float Q2);
  PhotonFlux(float W, float Q2, float beam_e);
  ~PhotonFlux();

  float GetVirtualPhotonFlux();
  std::string FluxString();
  void PrintFlux();
};

#endif
