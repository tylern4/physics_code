/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler												*/
/*	University Of South Carolina										*/
/************************************************************************/

#ifndef DELTA_T_CUT_H_GUARD
#define DELTA_T_CUT_H_GUARD
#include "main.h"
#include "delta_t_hist.hpp"

const double c_special_units = 29.9792458;
double vertex_time(double sc_time, double sc_pathlength, double cut_beta) {
	return sc_time - sc_pathlength/(cut_beta * c_special_units); 
}

double delta_t(double electron_vertex_time, double mass, double momentum, double sc_t, double sc_r){
	double cut_beta = 1.0/sqrt(1.0 + (mass/momentum)*(mass/momentum));
	return electron_vertex_time - vertex_time(sc_t,sc_r,cut_beta); 
}

void delta_t_cut(){
	double delta_t_P, delta_t_PIP, delta_t_ELECTRON;
	double electron_vertex = vertex_time(sc_t[sc[0]-1], sc_r[sc[0]-1], 1.0);

	for(int event_number = 0; event_number < gpart; event_number++){
		delta_t_P = delta_t(electron_vertex, MASS_P, p[event_number], sc_t[sc[event_number]-1], sc_r[sc[event_number]-1]);
		delta_t_PIP = delta_t(electron_vertex, MASS_PIP, p[event_number], sc_t[sc[event_number]-1], sc_r[sc[event_number]-1]);
		delta_t_ELECTRON = delta_t(electron_vertex, MASS_E, p[event_number], sc_t[sc[event_number]-1], sc_r[sc[event_number]-1]);

		delta_t_Fill(p[event_number], delta_t_P, delta_t_PIP, delta_t_ELECTRON);

		if((int)q[event_number] == 1) {
			Fill_deltat_P(p[event_number],delta_t_P);
			Fill_deltat_PIP(p[event_number],delta_t_PIP);
			Fill_deltat_positron(p[event_number], delta_t_ELECTRON);
			if(id[event_number] == PROTON) Fill_deltat_P_PID(p[event_number],delta_t_P);
			if(id[event_number] == PIP)    Fill_deltat_PIP_PID(p[event_number],delta_t_PIP);
			if(id[event_number] == -11)    Fill_deltat_positron_PID(p[event_number], delta_t_ELECTRON);
		} else if((int)q[event_number] == -1) {
			Fill_deltat_electron(p[event_number], delta_t_ELECTRON);
			if(id[event_number] == ELECTRON) Fill_deltat_electron_PID(p[event_number], delta_t_ELECTRON);
		}
	}
}
#endif
