/************************************************************************/
/*																		*/
/*																		*/
/*  Created by Nick Tyler												*/
/*	University Of South Carolina										*/
/************************************************************************/
#ifndef DELTA_TIME_H_GUARD
#define DELTA_TIME_H_GUARD
#include "delta_t.hpp"

class D_time
{
	double MASS_P = 0.93827203;
	double MASS_PIP = 0.13957018;
	double MASS_E = 0.000511;

	public:
		double proton_time = 0;
		double pip_time = 0;
		double electron_time = 0;
		double vertex_time = 0;
		double SC_T_vertex;
		double SC_R_vertex;
		double SC_T;
		double SC_R;
		double P;

	inline void SetVertexTimes(double sc_t_vertex, double sc_r_vertex){
		SC_T_vertex = sc_t_vertex;
		SC_R_vertex = sc_r_vertex;
	} //
	inline void SetTandP(double sc_t, double sc_r, double p){
		SC_T = sc_t;
		SC_R = sc_r;
		P = p;
	}

	inline D_time delta_t_calc(){
	
		Delta_T dt;
		D_time delta_time;

		double delta_t_P, delta_t_PIP, delta_t_ELECTRON;
		delta_time.vertex_time = dt.vertex_time(SC_T, SC_R, 1.0);

		delta_time.proton_time = dt.delta_t(delta_time.vertex_time, MASS_P, P, SC_T, SC_R);
		delta_time.pip_time = dt.delta_t(delta_time.vertex_time, MASS_PIP, P, SC_T, SC_R);
		delta_time.electron_time = dt.delta_t(delta_time.vertex_time, MASS_E, P, SC_T, SC_R);

		return delta_time;
	}


};
#endif