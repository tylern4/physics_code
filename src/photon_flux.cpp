/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#include "photon_flux.hpp"

PhotonFlux::PhotonFlux(float W, float Q2) : _W(W), _Q2(Q2) {
  _beam_momentum = Momentum(_beam_energy, MASS_E);
  _nu = photon_energy();
  _scattered_energy = (_beam_energy - _nu);
  _scattered_momentum = Momentum(_scattered_energy, MASS_E);

  _flux = photon_flux();
}

PhotonFlux::PhotonFlux(float W, float Q2, float beam_e) : _beam_energy(beam_e) { PhotonFlux(W, Q2); }

PhotonFlux::~PhotonFlux() {}

float PhotonFlux::Momentum(float E, float M) { return sqrtf(E * E - M * M); }

float PhotonFlux::photon_energy() { return ((_W * _W + _Q2) / _target_mass - _target_mass) / 2; }

float PhotonFlux::theta_calc() {
  return acosf((_beam_energy * _scattered_energy - _Q2 / 2.0 - MASS_E * MASS_E) /
               (_beam_momentum * _scattered_momentum));
}

float PhotonFlux::epsilon_calc() {
  return (powf((1 + (2 * (1 + ((_nu * _nu) / _Q2)) * powf(tanf((theta_calc() / 2)), 2))), -1));
}

float PhotonFlux::photon_flux() {
  return FS_ALPHA / (4 * PI * _Q2) * _W / (_beam_energy * _beam_energy * _target_mass * _target_mass) *
         (_W * _W - _target_mass * _target_mass) / (1 - epsilon_calc());
}

float PhotonFlux::GetVirtualPhotonFlux() { return _flux; }
std::string PhotonFlux::FluxString() { return std::to_string(_flux); }
void PhotonFlux::PrintFlux() { std::cout << RED << _flux << DEF << std::endl; }
