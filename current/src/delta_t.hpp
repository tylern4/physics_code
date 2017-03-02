/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef DELTA_T_H
#define DELTA_T_H
#include "histogram.hpp"
#include "constants.hpp"
#include "branches.hpp"

class Delta_T {
private:
  const double c_special_units = 29.9792458;

public:
  Delta_T();
  ~Delta_T();
  double vertex_time(double sc_time, double sc_pathlength, double cut_beta);
  double delta_t(double electron_vertex_time, double mass, double momentum,
                 double sc_t, double sc_r);
  void delta_t_cut(Histogram *hists, bool first_run);
  double *delta_t_array(double *dt_array, double mass);
  std::vector<double> delta_t_array(double mass, int num_parts);
};
#endif
