/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler												*/
/*	University Of South Carolina										*/
/************************************************************************/

#ifndef DELTA_T_H_GUARD
#define DELTA_T_H_GUARD

const double c_special_units = 29.9792458;

double vertex_time(double sc_time, double sc_pathlength, double cut_beta) {
	return sc_time - sc_pathlength/(cut_beta * c_special_units);
}

double delta_t(double electron_vertex_time, double mass, double momentum, double sc_t, double sc_r){
	double cut_beta = 1.0/sqrt(1.0 + pow(mass/momentum,2.0));
	return electron_vertex_time - vertex_time(sc_t,sc_r,cut_beta);
}
#endif