/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler												*/
/*	University Of South Carolina										*/
/************************************************************************/

#ifndef MISSING_H_GUARD
#define MISSING_H_GUARD
#include "TLorentzVector.h"

class MissingMass
{
	double MASS_P = 0.93827203;
	double MASS_PIP = 0.13957018;
	double MASS_E = 0.000511;

	double missing_mass_calc(TLorentzVector gamma_mu, TLorentzVector p_mu, TLorentzVector pi_mu){
		TLorentzVector reaction(0.0,0.0,0.0,0.0);
		reaction = (gamma_mu + p_mu - pi_mu);

		return reaction.M();
	} //

public:
	double mass = -1;
	double PX;
	double PY;
	double PZ;

	inline void MissingMassPxPyPz(double px, double py, double pz){
		PX = px;
		PY = py;
		PZ = pz;
	}

 	MissingMass missing_mass(TLorentzVector gamma_mu){
 		MissingMass MM;
		TVector3 Temp_vec_3_PIP(0.0,0.0,0.0);
		TLorentzVector Temp_vec_4_PIP(0.0,0.0,0.0,0.0);

		TLorentzVector Rest_PROTON(0.0,0.0,0.0,MASS_P);

		Temp_vec_3_PIP.SetXYZ(PX,PY,PZ);
		Temp_vec_4_PIP.SetVectM(Temp_vec_3_PIP, MASS_PIP);
		
		MM.mass = missing_mass_calc(gamma_mu,Rest_PROTON,Temp_vec_4_PIP);

		return MM;
	}

};
#endif
