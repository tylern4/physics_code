/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler												*/
/*	University Of South Carolina										*/
/************************************************************************/
#ifndef DELTA_T_H_GUARD
#define DELTA_T_H_GUARD
#include "TMath.h"

class Delta_T
{
	
public:
	inline double vertex_time(double sc_time, double sc_pathlength, double cut_beta) {
		return sc_time - sc_pathlength/(cut_beta * c_special_units); 
	} //

	inline double delta_t(double electron_vertex_time, double mass, double momentum, double sc_t, double sc_r){
		double cut_beta = 1.0/sqrt(1.0 + (mass/momentum)*(mass/momentum));
		return electron_vertex_time - vertex_time(sc_t,sc_r,cut_beta); 
	} //
private:
	const double c_special_units = 29.9792458;

};
#endif