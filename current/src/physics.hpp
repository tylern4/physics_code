/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#ifndef PHYSICS_H_GUARD
#define PHYSICS_H_GUARD
#include "constants.hpp"

namespace physics {
// Calcuating Q^2
// q^mu^2 = (e^mu - e^mu')^2 = -Q^2
double Q2_calc(LorentzVector e_mu, LorentzVector e_mu_prime);
//	Calcualting W
//	Gotten from s channel [(gamma - P)^2 == s == w^2]
//	Sqrt√[M_p^2 - Q^2 + 2 M_p gamma]
double W_calc(LorentzVector e_mu, LorentzVector e_mu_prime);
double xb_calc(double Q2, double E_prime);
double Q2_calc(const LorentzVector &gamma_mu);
double W_calc(const LorentzVector &gamma_mu);
// overload with 4 vectors instaed of other calculations
double xb_calc(LorentzVector e_mu, LorentzVector e_mu_prime);
double theta_calc(double cosz);
double phi_calc(double cosx, double cosy);
double center_phi_calc(double cosx, double cosy);
int get_sector(double phi);
double Get_Mass(int ID);
double fiducial_phi(double theta_e, double e_p);
LorentzVector *fourVec(double px, double py, double pz, double mass);
LorentzVector *fourVec(double p, double cx, double cy, double cz, double mass);
LorentzVector *fourVec(double px, double py, double pz, int pid);
LorentzVector *fourVec(double p, double cx, double cy, double cz, int pid);
float invTan(float y, float x);
float phi_boosted(std::unique_ptr<LorentzVector> &vec);

}  // namespace physics
#endif
