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
  float PX, PY, PZ;
  float target_mass, target_px, target_py, target_pz;
  float out_mass = MASS_PIP;
  float MM, MM2;
  TLorentzVector reaction;
  TLorentzVector vec_4_out;

 public:
  MissingMass();
  MissingMass(float t_mass, float t_p);
  ~MissingMass();

  void Set_PxPyPz(float px, float py, float pz);
  void Set_P_cos(float p_out, float cx_out, float cy_out, float cz_out);
  void Set_target_mass(float mass);
  void Set_target_PxPyPz(int zero);
  void Set_target_PxPyPz(float t_px, float t_py, float t_pz);
  void Set_target_P_cos(float t_p, float t_cx, float t_cy, float t_cz);
  void Set_4Vec(const TLorentzVector& event_p);
  void Add_4Vec(const TLorentzVector& event_p);
  void Reset();

  void missing_mass(TLorentzVector gamma_mu);
  float Get_MM();
  float Get_MM2();
};

#endif
