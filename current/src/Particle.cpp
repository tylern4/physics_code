/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "Particle.hpp"
// Constructor and destructor methods
Particle::Particle() {
  _vec4 = TLorentzVector(-99.0, -99.0, -99.0, -99.0);
  _vertex = TVector3(-99.0, -99.0, -99.0);
  _cos_angle = TVector3(-99.0, -99.0, -99.0);
}
Particle::Particle(double p, double cx, double cy, double cz, int pid) {
  _p = p;
  _PID = pid;
  _vertex = TVector3(-99.0, -99.0, -99.0);
  _mass = physics::Get_Mass(_PID);
  _cos_angle = TVector3(cx, cy, cz);
  _vec4 = physics::fourVec(_p, cx, cy, cz, _PID);
}
Particle::~Particle() {}

// Setup method in case you use the default constructer and fill later
void Particle::setup_4vector() { _vec4 = physics::fourVec(_p, _cos_angle.X(), _cos_angle.Y(), _cos_angle.Z(), _PID); }

// Set functions
void Particle::Set_momentum(double p, double cx, double cy, double cz) {
  _p = p;
  _cos_angle = TVector3(cx, cy, cz);
}
void Particle::Set_vertex(double vx, double vy, double vz) { _vertex = TVector3(vx, vy, vz); }
void Particle::Set_BankID(int pid) { _PID = pid; }
void Particle::Set_PID(int pid) { _bank_PID = pid; }
void Particle::Set_charge(int charge) { _charge = charge; }
void Particle::Set_beta(float beta) { _beta = beta; }
void Particle::Set_mass(float mass) { _mass = mass; }

// Get Function
TLorentzVector Particle::Particle4Vec() {
  if (_vec4.Px() != _vec4.Py() != _vec4.Pz() != _vec4.E()) {
    return _vec4;
  } else {
    setup_4vector();
    return _vec4;
  }
}

TVector3 Particle::vertex() { return _vertex; }

int Particle::PID() { return _PID; }
int Particle::BankID() { return _bank_PID; }
int Particle::charge() { return _charge; }
float Particle::beta() { return _beta; }
float Particle::mass() { return _mass; }
