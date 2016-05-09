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
//overload with 4 vectors instaed of otehr calculations
double xb_calc(TLorentzVector e_mu, TLorentzVector e_mu_prime){
	double Q2 = Q2_calc(e_mu,e_mu_prime);
	TLorentzVector q = e_mu - e_mu_prime;
	TLorentzVector target(0, 0, 0, MASS_P);
	return (Q2/ (2 * (q.Dot(target))));
}

//double missing_mass_calc(TLorentzVector gamma_mu, TLorentzVector p_mu, TLorentzVector pi_mu){
//	TLorentzVector reaction(0.0,0.0,0.0,0.0);
//	reaction = (gamma_mu + p_mu - pi_mu);
//
//	return reaction.M();
//}

double Get_Mass(int ID){

	switch (ID){
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
}

//	Print the readable name from particle ID
//	
void PrintID_Readable(int ID){
	switch (ID){
		case 2212:
			cout << "PROTON:";
			break;
		case 2112:
			cout << "NEUTRON:";
			break;
		case 211:
			cout << "PIP:";
			break;
		case -211:
			cout << "PIM:";
			break;
		case 111:
			cout << "PI0:";
			break;
		case 321:
			cout << "KP:";
			break;
		case -321:
			cout << "KM:";
			break;
		case 22:
			cout << "PHOTON:";
			break;
		case 11:
			cout << "ELECTRON:";
			break;
		case 0:
			cout << "***";
			break;
	}
}

//	Print the readable name from particle ID
//	
string PrintID_String(int ID){
	switch (ID){
		case 2212:
			return "PROTON";
			break;
		case 2112:
			return "NEUTRON";
			break;
		case 211:
			return "PIP";
			break;
		case -211:
			return "PIM";
			break;
		case 111:
			return "PI0";
			break;
		case 321:
			return "KP";
			break;
		case -321:
			return "KM";
			break;
		case 22:
			return "PHOTON";
			break;
		case 11:
			return "ELECTRON";
			break;
		case 0:
			return "***";
			break;
	}
}

#endif