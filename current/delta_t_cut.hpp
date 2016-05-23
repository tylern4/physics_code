/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler												*/
/*	University Of South Carolina										*/
/************************************************************************/

#ifndef DELTA_T_CUT_H_GUARD
#define DELTA_T_CUT_H_GUARD
#include "main.h"

const double c_special_units = 29.9792458;
double vertex_time(double sc_time, double sc_pathlength, double cut_beta) {
	return sc_time - sc_pathlength/(cut_beta * c_special_units); 
}

double delta_t(double electron_vertex_time, double mass, double momentum, double sc_t, double sc_r){
	double cut_beta = 1.0/sqrt(1.0 + (mass/momentum)*(mass/momentum));
	return electron_vertex_time - vertex_time(sc_t,sc_r,cut_beta); 
}

void delta_t_cut(){
	TVector3 Particle3(0.0,0.0,0.0);
	TLorentzVector Particle4(0.0,0.0,0.0,0.0);
	double delta_t_P, delta_t_PIP, delta_t_ELECTRON;
	double electron_vertex = vertex_time(sc_t[sc[0]-1], sc_r[sc[0]-1], 1.0);

	for(int event_number = 0; event_number < gpart; event_number++){
		//Get particles 3 and 4 vector for current event.
		Particle3.SetXYZ(p[event_number]*cx[event_number], p[event_number]*cy[event_number], p[event_number]*cz[event_number]);
		Particle4.SetVectM(Particle3,Get_Mass(id[event_number]));

		delta_t_P = delta_t(electron_vertex, MASS_P, p[event_number], sc_t[sc[event_number]-1], sc_r[sc[event_number]-1]);
		delta_t_PIP = delta_t(electron_vertex, MASS_PIP, p[event_number], sc_t[sc[event_number]-1], sc_r[sc[event_number]-1]);
		delta_t_ELECTRON = delta_t(electron_vertex, MASS_E, p[event_number], sc_t[sc[event_number]-1], sc_r[sc[event_number]-1]);

		if (event_number == 0 && id[0] == ELECTRON && gpart > 0 && stat[0] > 0 && (int)q[0] == -1 && sc[0] > 0 && dc[0] > 0 && ec[0] > 0 && dc_stat[dc[0]-1] > 0) {
			delta_t_Fill(p[event_number], delta_t_PIP, 8);
		}

		if (p[event_number] != 0 && (int)q[event_number] == 1) {
			delta_t_Fill(p[event_number],delta_t_P,3);

			delta_t_Fill(p[event_number], delta_t_PIP,4);

			delta_t_Fill(p[event_number], delta_t_ELECTRON, 5); 

			//If Pi+
			if(id[event_number] == PROTON && (int)q[event_number] == 1) {
				delta_t_Fill(p[event_number], delta_t_P, 1);
				delta_t_Fill(p[event_number], delta_t_ELECTRON, 7);
			//If Proton	
			} else if (id[event_number] == PIP && (int)q[event_number] == 1){
				delta_t_Fill(p[event_number], delta_t_PIP, 2);
				delta_t_Fill(p[event_number], delta_t_ELECTRON, 6);
			} 
		}
	} 

}
#endif
