/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#ifndef PHYSICS_H_GUARD
#define PHYSICS_H_GUARD
#include <cstdint>
#include <memory>
#include "Math/VectorUtil.h"
#include "constants.hpp"

namespace physics {
// Calcuating Q^2
// q^mu^2 = (e^mu - e^mu')^2 = -Q^2
double Q2_calc(const LorentzVector& e_mu, const LorentzVector& e_mu_prime);
//	Calcualting W
//	Gotten from s channel [(gamma - P)^2 == s == w^2]
//	Sqrt√[M_p^2 - Q^2 + 2 M_p gamma]
double W_calc(const LorentzVector& e_mu, const LorentzVector& e_mu_prime);
double xb_calc(double Q2, double E_prime);
double Q2_calc(const LorentzVector& gamma_mu);
double W_calc(const LorentzVector& gamma_mu);
// overload with 4 vectors instaed of other calculations
double xb_calc(const LorentzVector& gamma_mu);
double theta_rad(double cosz);
double phi_rad(double cosx, double cosy);
double theta_calc(double cosz);
double phi_calc(double cosx, double cosy);
double theta_calc_rad(double cosz);
double phi_calc_rad(double cosx, double cosy);
double center_phi_calc(double cosx, double cosy);
double center_phi_calc_rad(double cosx, double cosy);
int get_sector(double phi);
double Get_Mass(int ID);
double fiducial_phi(double theta_e, double e_p);
std::shared_ptr<LorentzVector> fourVec(double px, double py, double pz, double mass);
std::shared_ptr<LorentzVector> fourVec(double p, double cx, double cy, double cz, double mass);
std::shared_ptr<LorentzVector> fourVec(double px, double py, double pz, int pid);
std::shared_ptr<LorentzVector> fourVec(double p, double cx, double cy, double cz, int pid);
float invTan(const float& y, const float& x);
float phi_boosted(const std::shared_ptr<LorentzVector>& vec);

// True when `x` is a NaN or an infinity.
//
// These deliberately inspect the bit pattern instead of calling std::isnan /
// std::isinf. The build enables -Ofast (CMakeLists.txt), whose implied
// -ffinite-math-only lets the compiler assume no NaN or Inf is ever produced,
// so it is free to constant-fold those predicates to `false`. That silently
// disabled every NaN guard in the analysis: 51 of the 1524 output histograms
// were nondeterministic run to run, because NaN values were being binned into
// arbitrary bins (ROOT's bin selection for a NaN is unspecified), and a
// further 3 changed bin assignment on any recompile. A bit test is opaque to
// that optimisation, so these are always evaluated as written.
inline bool is_nan_or_inf(float x) {
  uint32_t bits;
  __builtin_memcpy(&bits, &x, sizeof(bits));
  return (bits & 0x7F800000u) == 0x7F800000u;  // exponent field all ones
}

inline bool is_finite(float x) { return !is_nan_or_inf(x); }

}  // namespace physics
#endif
