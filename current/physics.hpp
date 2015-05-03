/************************************************************************/
/*									
/*									
/*  Created by Nick Tyler					
/*	University Of South Carolina 😎			
/************************************************************************/

#ifndef PHYSICS_H_GUARD
#define PHYSICS_H_GUARD
#include "main.h"

using namespace std;

//	Calcualting W
//	Gotten from s channel [(gamma - P)^2 == s == w^2]
//	Sqrt√[M_p^2 - Q^2 + 2 M_p gamma]
double W_calc(double E_prime){
	//return sqrt( Square(MASS_P) + 2 * MASS_P * (E1D_E0-E_prime) );
	return sqrt(Square(MASS_P) - Q2 + 2 * MASS_P * (E1D_E0-E_prime));
}

//	Calulating Q^2 **Incorently
//	Gotten from t channel [(E_e - E_ep)^2 == t == -Q^2]
//	Q^2 = 4*E_beam*E_prime*Sin^2(theta/2)
double Q2_calc(double CosZ, double E_prime){
	double theta_2 = acos(CosZ)/2.0;
	double sin_sqr_theta_ovr_2 = Square(sin(theta_2));
	return 4 * E1D_E0 * E_prime * sin_sqr_theta_ovr_2;
}


//Calcuating Q^2 
// q^mu^2 = (e^mu - e^mu')^2 = -Q^2
double Q2_calc(TLorentzVector e_mu, TLorentzVector e_mu_prime){
	TLorentzVector q_mu = (e_mu - e_mu_prime);
	return q_mu.Mag2();
}

double xb_calc(double Q2, double E_prime){
	double gamma = E1D_E0-E_prime;
	double xb = (Q2/(2 * MASS_P * gamma));
	return xb;
}

double P_calc(double momentum, double CosX, double CosY, double CosZ){
	double Px = momentum * CosX;
	double Py = momentum * CosY;
	double Pz = momentum * CosZ;
	momentum = sqrt(Square(Px)+Square(Pz)+Square(Pz));
	return momentum;
}

//	Calcualting Energy from relativistic energy-momentum conservation
//	[E^2 = p^2 + m^2]
double E_calc(double momentum, double CosX, double CosY, double CosZ){
	momentum = P_calc(momentum,CosX,CosY,CosZ);
	double E2 = Square(momentum) + Square(MASS_E);

	return sqrt(E2);
}

//	Overloads of E_calc for different masses then electron
//	Calcualting Energy from relativistic energy-momentum conservation
//	[E^2 = p^2 + m^2]
double E_calc(double momentum, double CosX, double CosY, double CosZ, double mass){
	momentum = P_calc(momentum,CosX,CosY,CosZ);
	double E2 = Square(momentum) + Square(mass);

	return sqrt(E2);
}

//	Another overload with particle ID insead
double E_calc(double momentum, double CosX, double CosY, double CosZ, int ID){
	double mass;

	switch (ID){
		case 2212:
			mass = MASS_P;
			break;
		case 2112:
			mass = MASS_N;
			break;
		case 211:
			mass = MASS_PIP;
			break;
		case -211:
			mass = MASS_PIM;
			break;
		case 111:
			mass = MASS_PI0;
			break;
		case 321:
			mass = MASS_KP;
			break;
		case -321:
			mass = MASS_KM;
			break;
		case 22:
			mass = MASS_G;
			break;
		case 11:
			mass = MASS_E;
			break;
		case 0:
			mass = 0.0;
	}
	momentum = P_calc(momentum,CosX,CosY,CosZ);
	double E2 = Square(momentum)+Square(mass);

	return sqrt(E2);
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
	}
}

#endif