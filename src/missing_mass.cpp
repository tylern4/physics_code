/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "missing_mass.hpp"
#include <iostream>

MissingMass::MissingMass() {
  Set_target_mass(MASS_P);
  Set_target_PxPyPz(0);
  vec_4_out.SetXYZM(0, 0, 0, 0);
}
MissingMass::MissingMass(float t_mass, float t_p) {
  Set_target_mass(t_mass);
  Set_target_PxPyPz(t_p);
}
MissingMass::~MissingMass() {}

float MissingMass::Get_MM() { return MM; }
float MissingMass::Get_MM2() { return MM2; }

void MissingMass::Set_PxPyPz(float px, float py, float pz) {
  TVector3 vec_3_out;
  PX = px;
  PY = py;
  PZ = pz;
  vec_4_out.SetXYZM(PX, PY, PZ, out_mass);
}

void MissingMass::Set_4Vec(const TLorentzVector& event_p) { vec_4_out = event_p; }
void MissingMass::Add_4Vec(const TLorentzVector& event_p) { vec_4_out += event_p; }
void MissingMass::Reset() { vec_4_out.SetXYZM(0, 0, 0, 0); }

void MissingMass::Set_P_cos(float p_out, float cx_out, float cy_out, float cz_out) {
  PX = p_out * cx_out;
  PY = p_out * cy_out;
  PZ = p_out * cz_out;
}

void MissingMass::Set_target_mass(float mass) { target_mass = mass; }

void MissingMass::Set_target_PxPyPz(int zero) {
  target_px = 0.0;
  target_py = 0.0;
  target_pz = 0.0;
}

void MissingMass::Set_target_PxPyPz(float t_px, float t_py, float t_pz) {
  target_px = t_px;
  target_py = t_py;
  target_pz = t_pz;
}

void MissingMass::Set_target_P_cos(float t_p, float t_cx, float t_cy, float t_cz) {
  target_px = t_p * t_cx;
  target_py = t_p * t_cy;
  target_pz = t_p * t_cz;
}

void MissingMass::missing_mass(TLorentzVector gamma_mu) {
  // Initialize all vectors
  TVector3 target_3;
  TLorentzVector target;

  // Set target vector
  target_3.SetXYZ(target_px, target_py, target_pz);
  target.SetVectM(target_3, target_mass);
  // Set output particle vector

  reaction = (gamma_mu + target - vec_4_out);

  MM2 = reaction.M2();
  MM = reaction.M();
}
