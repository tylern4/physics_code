/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler												*/
/*	University Of South Carolina										*/
/************************************************************************/

#ifndef PHYSICS_H_GUARD
#define PHYSICS_H_GUARD
#include "main.h"
#include <TLorentzVector.h>
#include "TROOT.h"

//Calcuating Q^2
// q^mu^2 = (e^mu - e^mu')^2 = -Q^2
double Q2_calc(TLorentzVector e_mu, TLorentzVector e_mu_prime){
	TLorentzVector q_mu = (e_mu - e_mu_prime);
	return -q_mu.Mag2();
}
//	Calcualting W
//	Gotten from s channel [(gamma - P)^2 == s == w^2]
//	Sqrt√[M_p^2 - Q^2 + 2 M_p gamma]
double W_calc(TLorentzVector e_mu, TLorentzVector e_mu_prime){
	TLorentzVector q_mu = (e_mu - e_mu_prime);
	TVector3 p_mu_3(0,0,0);
	TLorentzVector p_mu;
	p_mu.SetVectM(p_mu_3,MASS_P);
	return (p_mu + q_mu).Mag();
}

double xb_calc(double Q2, double E_prime){
	double gamma = E1D_E0-E_prime;
	double xb = (Q2/(2 * MASS_P * gamma));
	return xb;
}
//overload with 4 vectors instaed of other calculations
double xb_calc(TLorentzVector e_mu, TLorentzVector e_mu_prime){
	double Q2 = Q2_calc(e_mu,e_mu_prime);
	TLorentzVector q = e_mu - e_mu_prime;
	TLorentzVector target(0, 0, 0, MASS_P);
	return (Q2/ (2 * (q.Dot(target))));
}

double theta_calc(double cosz){
	return acos(cosz)/D2R;
}

double phi_calc(double cosx, double cosy){
	return atan2(cosx, cosy)/D2R;
}

double center_phi_calc(double cosx, double cosy){
	double phi0 = (atan2(cosx, cosy)/D2R);
	phi0 += 30;
 	if (phi0 < 0.0) phi0 += 360.0;
	if (phi0 > 360.0) phi0 -= 360.0;
	return phi0;
}

int get_sector(double phi) {
	if(phi>=0 && phi <60) {
		return 0;
	} else if(phi>=60 && phi<120) {
		return 1;
	} else if(phi>=120 && phi <=180) {
		return 2;
	} else if(phi>=-180 && phi<-120) {
		return 3;
	} else if(phi>=-120 && phi<-60) {
		return 4;
	} else if(phi>=-60 && phi<0) {
		return 5;
	} else {
		return (int)std::nan("0");
	}
}


#endif
