/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#include "delta_t.hpp"

Delta_T::Delta_T() {}
Delta_T::Delta_T(double sc_time, double sc_pathlength) { vertex = vertex_time(sc_time, sc_pathlength, 1.0); }

Delta_T::~Delta_T() {}

double Delta_T::vertex_time(double sc_time, double sc_pathlength, double relatavistic_beta) {
  return sc_time - sc_pathlength / (relatavistic_beta * c_special_units);
}

void Delta_T::deltat(double momentum, double sc_t, double sc_r) {
  double beta = 0.0;
  double mp = (masses.at(0) / momentum);
  beta = 1.0 / sqrt(1.0 + (mp * mp));
  dt_E = vertex - vertex_time(sc_t, sc_r, beta);

  mp = (masses.at(1) / momentum);
  beta = 1.0 / sqrt(1.0 + (mp * mp));
  dt_P = vertex - vertex_time(sc_t, sc_r, beta);

  mp = (masses.at(2) / momentum);
  beta = 1.0 / sqrt(1.0 + (mp * mp));
  dt_Pi = vertex - vertex_time(sc_t, sc_r, beta);

  mp = (masses.at(3) / momentum);
  beta = 1.0 / sqrt(1.0 + (mp * mp));
  dt_K = vertex - vertex_time(sc_t, sc_r, beta);
}

double Delta_T::deltat(double momentum, double sc_t, double sc_r, double mass) {
  double beta = 0.0;
  double mp = (mass / momentum);
  beta = 1.0 / sqrt(1.0 + (mp * mp));
  double out = vertex - vertex_time(sc_t, sc_r, beta);
  return out;
}

double Delta_T::Get_dt_E() { return dt_E; }
double Delta_T::Get_dt_P() { return dt_P; }
double Delta_T::Get_dt_Pi() { return dt_Pi; }
double Delta_T::Get_dt_K() { return dt_K; }
double Delta_T::Get_vertex() { return vertex; }

double Delta_T::delta_t(double electron_vertex_time, double mass, double momentum, double sc_t, double sc_r) {
  double cut_beta = 1.0 / sqrt(1.0 + (mass / momentum) * (mass / momentum));
  return electron_vertex_time - vertex_time(sc_t, sc_r, cut_beta);
}
