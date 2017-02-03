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

double *delta_t_array(double* dt_array, double mass){
	double delta_t_P;
	double electron_vertex = vertex_time(sc_t[sc[0]-1], sc_r[sc[0]-1], 1.0);
	double sct,scr,mom;
	int ID,charge,sc_paddle,sc_sector;

	for(int event_number = 0; event_number < gpart; event_number++){
		sct = (double)sc_t[sc[event_number]-1];
		scr = (double)sc_r[sc[event_number]-1];
		mom = (double)p[event_number];
		ID = (int)id[event_number];
		charge = (int)q[event_number];
		sc_paddle = (int)sc_pd[sc[event_number]-1];
		sc_sector = (int)sc_sect[sc[event_number]-1];

		dt_array[event_number] = delta_t(electron_vertex, mass, mom, sct, scr);
	}
	return dt_array;
}

std::vector<double> delta_t_array(double mass, int num_parts){
	std::vector<double> dt_array(num_parts);
	double delta_t_P;
	double electron_vertex = vertex_time(sc_t[sc[0]-1], sc_r[sc[0]-1], 1.0);
	double sct,scr,mom;
	int ID,charge,sc_paddle,sc_sector;

	for(int event_number = 0; event_number < gpart; event_number++){
		sct = (double)sc_t[sc[event_number]-1];
		scr = (double)sc_r[sc[event_number]-1];
		mom = (double)p[event_number];
		ID = (int)id[event_number];
		charge = (int)q[event_number];
		sc_paddle = (int)sc_pd[sc[event_number]-1];
		sc_sector = (int)sc_sect[sc[event_number]-1];

		dt_array[event_number] = delta_t(electron_vertex, mass, mom, sct, scr);
	}
	return dt_array;
}

#endif
