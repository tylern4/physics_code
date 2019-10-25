/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "physics.hpp"

std::shared_ptr<LorentzVector> physics::fourVec(double px, double py, double pz, double mass) {
  return std::make_shared<LorentzVector>(px, py, pz, mass);
}
std::shared_ptr<LorentzVector> physics::fourVec(double px, double py, double pz, int pid) {
  return physics::fourVec(px, py, pz, physics::Get_Mass(pid));
}
std::shared_ptr<LorentzVector> physics::fourVec(double p, double cx, double cy, double cz, double mass) {
  return physics::fourVec(p * cx, p * cy, p * cz, mass);
}

std::shared_ptr<LorentzVector> physics::fourVec(double p, double cx, double cy, double cz, int pid) {
  return physics::fourVec(p * cx, p * cy, p * cz, pid);
}

// Calcuating Q^2
//	Gotten from t channel
// -q^mu^2 = -(e^mu - e^mu')^2 = Q^2
double physics::Q2_calc(const LorentzVector &e_mu, const LorentzVector &e_mu_prime) {
  LorentzVector q_mu = (e_mu - e_mu_prime);
  return -q_mu.mag2();
}
//	Calcualting W
//	Gotten from s channel [(gamma + P)^2 == s == w^2]
//	Sqrt√[M_p^2 - Q^2 + 2 M_p gamma]
double physics::W_calc(const LorentzVector &e_mu, const LorentzVector &e_mu_prime) {
  LorentzVector q_mu = (e_mu - e_mu_prime);
  LorentzVector p_mu(0, 0, 0, MASS_P);
  return (p_mu + q_mu).mag();
}

double physics::Q2_calc(const LorentzVector &gamma_mu) { return -gamma_mu.mag2(); }

double physics::W_calc(const LorentzVector &gamma_mu) {
  LorentzVector p_mu(0.0, 0.0, 0.0, MASS_P);
  return (p_mu + gamma_mu).mag();
}

double physics::xb_calc(double Q2, double E_prime) {
  double gamma = E1D_E0 - E_prime;
  double xb = (Q2 / (2 * MASS_P * gamma));
  return xb;
}

double physics::xb_calc(const LorentzVector &gamma_mu) {
  double Q2 = Q2_calc(gamma_mu);
  LorentzVector target(0, 0, 0, MASS_P);
  return (Q2 / (2 * (gamma_mu.Dot(target))));
}

double physics::theta_calc(double cosz) { return acos(cosz) * RAD2DEG; }

double physics::phi_calc(double cosx, double cosy) { return atan2(cosx, cosy) * RAD2DEG; }

float physics::invTan(float y, float x) {
  if (x > 0 && y > 0)
    return atan(y / x);  // 1st Quad.
  else if (x < 0 && y > 0)
    return atan(y / x) + PI;  // 2nd Quad
  else if (x < 0 && y < 0)
    return atan(y / x) + PI;  // 3rd Quad
  else if (x > 0 && y < 0)
    return atan(y / x) + 2 * PI;  // 4th Quad
  else if (x == 0 && y > 0)
    return PI / 2;
  else if (x == 0 && y < 0)
    return 3 * PI / 2;
  return NAN;
}

float physics::phi_boosted(const std::shared_ptr<LorentzVector> &vec) { return invTan(vec->Py(), vec->Px()); }

double physics::center_phi_calc(double cosx, double cosy) {
  double phi0 = (atan2(cosx, cosy) * RAD2DEG);
  phi0 += 30;
  if (phi0 < 0.0) phi0 += 360.0;
  if (phi0 > 360.0) phi0 -= 360.0;
  return phi0;
}

int physics::get_sector(double phi) {
  if (phi >= 60 && phi < 120) {
    return 1;
  } else if (phi >= 0 && phi < 60) {
    return 2;
  } else if (phi >= -60 && phi < 0) {
    return 3;
  } else if (phi >= -120 && phi < -60) {
    return 4;
  } else if (phi >= -180 && phi < -120) {
    return 5;
  } else if (phi >= 120 && phi < 180) {
    return 6;
  } else {
    return (int)NULL;
  }
}

double physics::Get_Mass(int ID) {
  switch (ID) {
    case 2212:
      return MASS_P;
      break;
    case 2112:
      return MASS_N;
      break;
    case 211:
      return MASS_PIP;
      break;
    case -211:
      return MASS_PIM;
      break;
    case 111:
      return MASS_PI0;
      break;
    case 321:
      return MASS_KP;
      break;
    case -321:
      return MASS_KM;
      break;
    case 22:
      return MASS_G;
      break;
    case 11:
      return MASS_E;
      break;
    case 0:
      return 0.0;
      break;
  }
  return 0;
}
