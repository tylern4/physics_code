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

using namespace std;

//Calcuating Q^2 
// q^mu^2 = (e^mu - e^mu')^2 = -Q^2
double Q2_calc(TLorentzVector e_mu, TLorentzVector e_mu_prime){
	TLorentzVector q_mu = (e_mu - e_mu_prime);
	return -q_mu.Mag2();
}

//	Calcualting W
//	Gotten from s channel [(gamma - P)^2 == s == w^2]
//	Sqrt√[M_p^2 - Q^2 + 2 M_p gamma]
/*double W_calc(double E_prime){
	//return sqrt( Square(MASS_P) + 2 * MASS_P * (E1D_E0-E_prime) );
	return sqrt(Square(MASS_P) - Q2 + 2 * MASS_P * (E1D_E0-E_prime));
}*/

/*double W_calc(TLorentzVector e_mu, TLorentzVector e_mu_prime){
	return sqrt(Square(MASS_P) - Q2_calc(e_mu, e_mu_prime) + 2 * MASS_P * (e_mu.E() - e_mu_prime.E()));
}*/

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

/*double P_calc(double momentum, double CosX, double CosY, double CosZ){
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
*/
//	Another overload with particle ID insead
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
	}
}

#endif