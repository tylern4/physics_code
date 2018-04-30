/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "particle.hpp"

Particle::Particle() { vec = new TLorentzVector; }
Particle::Particle(double p, double cx, double cy, double cz, int id, double m)
    : px(p * cx), py(p * cy), pz(p * cz), pid(id), mass(m) {
  vec = new TLorentzVector;
  theta = physics::theta_calc(cz);
  phi = physics::phi_calc(cx, cy);
  sector = physics::get_sector(phi);
  *vec = physics::fourVec(px, py, pz, mass);
}
Particle::Particle(double p, double cx, double cy, double cz, int id) : px(p * cx), py(p * cy), pz(p * cz), pid(id) {
  vec = new TLorentzVector;
  mass = physics::Get_Mass(pid);
  theta = physics::theta_calc(cz);
  phi = physics::phi_calc(cx, cy);
  sector = physics::get_sector(phi);
  *vec = physics::fourVec(px, py, pz, physics::Get_Mass(pid));
}
Particle::~Particle() { delete vec; }

inline void Particle::ReCalc() {
  theta = physics::theta_calc(cz);
  phi = physics::phi_calc(cx, cy);
  sector = physics::get_sector(phi);
  *vec = physics::fourVec(px, py, pz, mass);
}

inline double Particle::SetPx(double _px) { px = _px; }
inline double Particle::SetPy(double _py) { py = _py; }
inline double Particle::SetPz(double _pz) { pz = _pz; }
inline double Particle::SetP(double _p) { p = _p; }
inline double Particle::SetMass(double _mass) { mass = _mass; }
inline double Particle::SetTheta(double _theta) { theta = _theta; }
inline double Particle::SetPhi(double _phi) { phi = _phi; }
inline int Particle::SetSector(int _sec) { sector = _sec; }
inline int Particle::SetPid(int _pid) { pid = _pid; }

double Particle::GetPx() { return px; }
double Particle::GetPy() { return py; }
double Particle::GetPz() { return pz; }
double Particle::GetP() { return p; }
double Particle::GetE() { return vec->E(); }
double Particle::GetMass() { return mass; }
double Particle::GetTheta() { return theta; }
double Particle::GetPhi() { return phi; }
int Particle::GetSector() { return sector; }
int Particle::GetPid() { return pid; }
TLorentzVector Particle::Get4Vec() { return *vec; }

Electron::Electron(double p, double cx, double cy, double cz) : Particle(p, cx, cy, cz, ELECTRON, MASS_E) {
  *q_mu = *beam - this->Get4Vec();
  W = physics::W_calc(this->Get4Vec());
  Q2 = physics::Q2_calc(this->Get4Vec());
}
Electron::~Electron() {
  delete beam;
  delete q_mu;
}

double Electron::GetW() { return W; }
double Electron::GetQ2() { return Q2; }
TLorentzVector Electron::Get_q_mu() { return *q_mu; }
