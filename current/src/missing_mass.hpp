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
  double PX, PY, PZ;
  double target_mass, target_px, target_py, target_pz;
  double out_mass = MASS_PIP;
  double MM, MM2;
  TLorentzVector reaction;
  TLorentzVector vec_4_out;

 public:
  MissingMass();
  MissingMass(double t_mass, double t_p);
  ~MissingMass();

  void Set_PxPyPz(double px, double py, double pz);
  void Set_P_cos(double p_out, double cx_out, double cy_out, double cz_out);
  void Set_target_mass(double mass);
  void Set_target_PxPyPz(int zero);
  void Set_target_PxPyPz(double t_px, double t_py, double t_pz);
  void Set_target_P_cos(double t_p, double t_cx, double t_cy, double t_cz);
  void Set_4Vec(TLorentzVector event_p);
  void Add_4Vec(TLorentzVector event_p);
  void Reset();

  void missing_mass(TLorentzVector gamma_mu);
  double Get_MM();
  double Get_MM2();
};

#endif
