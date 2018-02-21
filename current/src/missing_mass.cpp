/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "missing_mass.hpp"

MissingMass::MissingMass() {
  Set_target_mass(MASS_P);
  Set_target_PxPyPz(0);
}
MissingMass::MissingMass(double t_mass, double t_p) {
  Set_target_mass(t_mass);
  Set_target_PxPyPz(t_p);
}
MissingMass::~MissingMass() {}

double MissingMass::Get_MM() { return MM; }
double MissingMass::Get_MM2() { return MM2; }

void MissingMass::Set_PxPyPz(double px, double py, double pz) {
  PX = px;
  PY = py;
  PZ = pz;
}

void MissingMass::Set_P_cos(double p_out, double cx_out, double cy_out, double cz_out) {
  PX = p_out * cx_out;
  PY = p_out * cy_out;
  PZ = p_out * cz_out;
}

void MissingMass::Set_target_mass(double mass) { target_mass = mass; }

void MissingMass::Set_target_PxPyPz(int zero) {
  target_px = 0.0;
  target_py = 0.0;
  target_pz = 0.0;
}

void MissingMass::Set_target_PxPyPz(double t_px, double t_py, double t_pz) {
  target_px = t_px;
  target_py = t_py;
  target_pz = t_pz;
}

void MissingMass::Set_target_P_cos(double t_p, double t_cx, double t_cy, double t_cz) {
  target_px = t_p * t_cx;
  target_py = t_p * t_cy;
  target_pz = t_p * t_cz;
}

void MissingMass::missing_mass(TLorentzVector gamma_mu) {
  // Initialize all vectors
  TVector3 vec_3_out;
  TLorentzVector vec_4_out;
  TVector3 target_3;
  TLorentzVector target;

  // Set target vector
  target_3.SetXYZ(target_px, target_py, target_pz);
  target.SetVectM(target_3, target_mass);
  // Set output particle vector
  vec_3_out.SetXYZ(PX, PY, PZ);
  vec_4_out.SetVectM(vec_3_out, out_mass);

  reaction = (gamma_mu + target - vec_4_out);

  MM2 = reaction.M2();
  MM = reaction.M();
}
