/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#ifndef MISSING_H
#define MISSING_H
#include "TLorentzVector.h"
#include "constants.hpp"

class MissingMass {
 private:
  double PX;
  double PY;
  double PZ;
  double target_mass = MASS_P;
  double target_px = 0.0;
  double target_py = 0.0;
  double target_pz = 0.0;
  double out_mass = MASS_PIP;

 public:
  MissingMass();
  ~MissingMass();

  void Set_PxPyPz(double px, double py, double pz);
  void Set_P_cos(double p_out, double cx_out, double cy_out, double cz_out);
  void Set_target_mass(double mass);
  void Set_target_PxPyPz(int zero);
  void Set_target_PxPyPz(double t_px, double t_py, double t_pz);
  void Set_target_P_cos(double t_p, double t_cx, double t_cy, double t_cz);

  double missing_mass(TLorentzVector gamma_mu);
};

#endif
