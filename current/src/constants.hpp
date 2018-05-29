/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#ifndef CONSTANTS_H_GUARD
#define CONSTANTS_H_GUARD
#include "TMath.h"

static const int MAX_PARTS = 100;
static const int N_SIGMA = 3;
static const double PI = TMath::Pi();
static const float INV_SQRT_2PI = TMath::Sqrt(2 * TMath::Pi());
static const double D2R = PI / 180.0;
// static const double E1D_E0 = 4.802;   // GeV
// Run Info 22848 - 23500
// http://clasweb.jlab.org/clasonline/servlet/runinfo?action=detail&run=22880
static const double E1D_E0 = 4.81726;  // GeV

static const double SOL = 29.9792458;
// misc. constants
static const double FSC = 0.00729735253;
static const double NA = 6.02214129E23;   // Avigadro's number
static const double QE = 1.60217646E-19;  // Charge or electron

// particle codes, usually PDG codes, but always those used in BOS
static const int PROTON = 2212;
static const int NEUTRON = 2112;
static const int PIP = 211;
static const int PIM = -211;
static const int PI0 = 111;
static const int KP = 321;
static const int KM = -321;
static const int PHOTON = 22;
static const int ELECTRON = 11;

// PDG particle masses in GeV/c2
static const double MASS_P = 0.93827203;
static const double MASS_N = 0.93956556;
static const double MASS_E = 0.000511;
static const double MASS_PIP = 0.13957018;
static const double MASS_PIM = 0.13957018;
static const double MASS_PI0 = 0.1349766;
static const double MASS_KP = 0.493677;
static const double MASS_KM = 0.493677;
static const double MASS_G = 0.0;
static const double MASS_OMEGA = 0.78265;

// Got from ye/arjun
// TODO::Get right EC threshhold
// p_min(in MeV) = 214 + 2.47 × EC_threshold(in mV)
static const double MIN_P_CUT = 0.7000;

const double min_phi[6] = {0, 60, 120, -180, -120, -60};
const double max_phi[6] = {60, 120, 180, -120, -60, 0};

#endif
