/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/
#include "physics.hpp"

// Calcuating Q^2
// q^mu^2 = (e^mu - e^mu')^2 = -Q^2
double Q2_calc(TLorentzVector e_mu, TLorentzVector e_mu_prime) {
  TLorentzVector q_mu = (e_mu - e_mu_prime);
  return -q_mu.Mag2();
}
//	Calcualting W
//	Gotten from s channel [(gamma - P)^2 == s == w^2]
//	Sqrt√[M_p^2 - Q^2 + 2 M_p gamma]
double W_calc(TLorentzVector e_mu, TLorentzVector e_mu_prime) {
  TLorentzVector q_mu = (e_mu - e_mu_prime);
  TVector3 p_mu_3(0, 0, 0);
  TLorentzVector p_mu;
  p_mu.SetVectM(p_mu_3, MASS_P);
  return (p_mu + q_mu).Mag();
}

double xb_calc(double Q2, double E_prime) {
  double gamma = E1D_E0 - E_prime;
  double xb = (Q2 / (2 * MASS_P * gamma));
  return xb;
}
// overload with 4 vectors instaed of otehr calculations
double xb_calc(TLorentzVector e_mu, TLorentzVector e_mu_prime) {
  double Q2 = Q2_calc(e_mu, e_mu_prime);
  TLorentzVector q = e_mu - e_mu_prime;
  TLorentzVector target(0, 0, 0, MASS_P);
  return (Q2 / (2 * (q.Dot(target))));
}

double theta_calc(double cosz) { return acos(cosz) / D2R; }

double phi_calc(double cosx, double cosy) { return atan2(cosx, cosy) / D2R; }

double center_phi_calc(double cosx, double cosy) {
  double phi0 = (atan2(cosx, cosy) / D2R);
  phi0 += 30;
  if (phi0 < 0.0)
    phi0 += 360.0;
  if (phi0 > 360.0)
    phi0 -= 360.0;
  return phi0;
}

int get_sector(double phi) {
  if (phi >= 0 && phi < 60) {
    return 0;
  } else if (phi >= 60 && phi < 120) {
    return 1;
  } else if (phi >= 120 && phi <= 180) {
    return 2;
  } else if (phi >= -180 && phi < -120) {
    return 3;
  } else if (phi >= -120 && phi < -60) {
    return 4;
  } else if (phi >= -60 && phi < 0) {
    return 5;
  } else {
    return (int)std::nan("0");
  }
  /*if(phi>=-30 && phi <30) {
          return 0;
  } else if(phi>=30 && phi<90) {
          return 1;
  } else if(phi>=90 && phi <150) {
          return 2;
  } else if(phi>=150 || phi<-150) {
          return 3;
  } else if(phi>=-150 && phi<-90) {
          return 4;
  } else if(phi>=-90 && phi<-30) {
          return 5;
  } else {
          return (int)std::nan("0");
  } */
}

double Get_Mass(int ID) {

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