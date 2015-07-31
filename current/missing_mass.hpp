/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler												*/
/*	University Of South Carolina										*/
/************************************************************************/

#ifndef MISSING_H_GUARD
#define MISSING_H_GUARD
#include <stdlib.h>
#include <stdio.h>
#include "main.h"

void missing_mass(TLorentzVector gamma_mu){
	TVector3 Temp_vec_3_PIP(0.0,0.0,0.0);
	TLorentzVector Temp_vec_4_PIP(0.0,0.0,0.0,0.0);

	TVector3 Temp_vec_3_PROTON(0.0,0.0,0.0);
	TLorentzVector Temp_vec_4_PROTON(0.0,0.0,0.0,MASS_P);
	int numOfPis = 0, numOfProtons = 0;

	for(int event_number = 0; event_number < gpart; event_number++) {
		/*if (id[event_number] == PROTON && (int)q[event_number] == 1) {
			Temp_vec_3_PROTON.SetXYZ(p[event_number]*cx[event_number], p[event_number]*cy[event_number], p[event_number]*cz[event_number]);
			Temp_vec_4_PROTON.SetVectM(Temp_vec_3_PROTON, MASS_P);
			numOfProtons++;
		} else*/ if (id[event_number] == PIP && (int)q[event_number] == 1) {
			numOfPis++;
			Temp_vec_3_PIP.SetXYZ(p[event_number]*cx[event_number], p[event_number]*cy[event_number], p[event_number]*cz[event_number]);
			Temp_vec_4_PIP.SetVectM(Temp_vec_3_PIP, MASS_PIP);
		} 
	}
	if (numOfPis == 1 /*&& numOfProtons == 1*/) Fill_Missing_Mass(missing_mass_calc(gamma_mu,Temp_vec_4_PROTON,Temp_vec_4_PIP));
	if (numOfPis == 1 && numOfProtons == 1 && gpart == 3) {
		/*cout << red << "************" << endl;
		cout << "gamma: ";
		gamma_mu.Print();
		cout << "proton: ";
		Temp_vec_4_PROTON.Print();
		cout << "pi: ";
		Temp_vec_4_PIP.Print();
		cout << green << missing_mass_calc(gamma_mu,Temp_vec_4_PROTON,Temp_vec_4_PIP) << endl; */

		Fill_Missing_Mass_P_PI(missing_mass_calc(gamma_mu,Temp_vec_4_PROTON,Temp_vec_4_PIP));
	}
}

#endif
