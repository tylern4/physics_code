/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#ifndef CONSTANTS_H_GUARD
#define CONSTANTS_H_GUARD
#include <iostream>
#include <unordered_map>
#include "Math/Rotation3D.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TMath.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<float>> LorentzVector;
typedef ROOT::Math::Rotation3D LorentzRotation;
typedef ROOT::Math::XYZVectorD Vector3;

static const short NUM_THREADS = (getenv("NUM_THREADS") != NULL) ? atoi(getenv("NUM_THREADS")) : 8;
// static const float BEAM_E = (getenv("BEAM_E") != NULL) ? atof(getenv("BEAM_E")) : 4.81726;

static const short NUM_SECTORS = 6;
static const short MAX_PARTS = 100;
static const short N_SIGMA = 3;
static const double PI = TMath::Pi();
static const double INV_SQRT_2PI = TMath::Sqrt(2 * TMath::Pi());
static const double D2R = PI / 180.0;
static const double DEG2RAD = PI / 180.0;
static const double RAD2DEG = 180.0 / PI;
static const int POSITIVE = 1;
static const int NEGATIVE = -1;
// static const double E1D_E0 = 4.802;   // GeV
// Run Info 22848 - 23500
// http://clasweb.jlab.org/clasonline/servlet/runinfo?action=detail&run=22880
static const double E1D_E0 = 4.81726;  // GeV
static const double E1F_E0 = 5.479;    // GeV

static const double SOL = 29.9792458;
// misc. constants
static const double FSC = 0.00729735253;
static const double NA = 6.02214129E23;               // Avigadro's number
static const double QE = 1.60217646E-19;              // Charge or electron
static const double FS_ALPHA = 0.007297352570866302;  // Fine structure alpha

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
static const double MIN_P_CUT = 0.005000;

// const double min_phi[6] = {0, 60, 120, -180, -120, -60};
// const double max_phi[6] = {60, 120, 180, -120, -60, 0};

static std::unordered_map<int, std::pair<float, float>> phi_map = {
    {1, std::make_pair(60, 120)},   {2, std::make_pair(0, 60)},      {3, std::make_pair(-60, 0)},
    {4, std::make_pair(-120, -60)}, {5, std::make_pair(-180, -120)}, {6, std::make_pair(120, 180)}};

static std::unordered_map<int, float> phi_center = {{1, 90}, {2, 30}, {3, -30}, {4, -90}, {5, -150}, {6, 150}};

const double min_phi[6] = {60, 0, -60, -120, -180, 120};
const double max_phi[6] = {120, 60, 0, -60, -120, 180};

static std::unordered_map<int, double> mass_map = {
    {PROTON, MASS_P}, {-PROTON, MASS_P}, {NEUTRON, MASS_N}, {PIP, MASS_PIP},    {PIM, MASS_PIM},    {PI0, MASS_PI0},
    {KP, MASS_KP},    {KM, MASS_KM},     {PHOTON, MASS_G},  {ELECTRON, MASS_E}, {-ELECTRON, MASS_E}};

static std::unordered_map<int, std::string> particle_name = {
    {PROTON, "P"}, {-PROTON, "-P"}, {NEUTRON, "N"}, {PIP, "PIP"},    {PIM, "PIM"},    {PI0, "PI0"},
    {KP, "KP"},    {KM, "KM"},      {PHOTON, "G"},  {ELECTRON, "E"}, {-ELECTRON, "E"}};

static const float a0xh[6] = {24.0, 24.0, 23.0, 23.5, 24.5, 24.5};
static const float a0mh[6] = {25.0, 26.0, 26.0, 25.5, 27.0, 26.0};
static const float a1xh[6] = {0.22, 0.23, 0.20, 0.20, 0.22, 0.22};
static const float a1mh[6] = {0.22, 0.22, 0.22, 0.22, 0.16, 0.16};
static const float a2xh[6] = {8.0, 8.0, 8.0, 8.0, 8.0, 8.0};
static const float a2mh[6] = {8.0, 8.0, 8.0, 8.0, 8.0, 8.0};
static const float a3xh[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
static const float a3mh[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

static const std::vector<double> dt_pip_const_top = {0.026476496274405754, -0.20219089310509566, 0.5065924067561807,
                                                     -0.4949841169578805, 0.5621128407552162};
static const std::vector<double> dt_pip_const_bottom = {-0.028629453054453677, 0.23665267204580884, -0.6637868377197707,
                                                        0.7411629911184183, -0.6968506895457166};
static const std::vector<double> dt_P_const_top = {0.08240030048909537, -0.7710826520788646, 2.570661547673607,
                                                   -3.509283013250415, 1.9410227467899008};
static const std::vector<double> dt_P_const_bottom = {-0.1433929584211291, 1.2741941832017054, -4.027205118341501,
                                                      5.306509383159595, -2.863560530067931};

static const std::vector<double> dt_pip_top = {-0.5680719469021316, 4.065155482173403,  -2.285775227566277,
                                               0.09291303310648173, 1.6698856065697703, 0.6181398119277808};

static const std::vector<double> dt_pip_bottom = {0.7021097407492168, 1.1754282714198145,  0.851129051933132,
                                                  0.8804059391935983, -1.9632091423474742, 0.22761907815192178};

enum mm2_names { A, mu, sigma, lambda1, lambda2 };
static const std::vector<std::vector<double>> mm2_fit_values = {
    {40.04364142911593660301, 0.89586494546690320639, 0.01421356258481163565, 54.12038299108210281929,
     19.88478985124503850557},
    {24.89731359731760917953, 0.87883493195336415127, 0.01302611651420351735, 48.53601194230342485980,
     26.91075411204556644407},
    {26.92573964785830753499, 0.90531823888385554167, 0.01317927668781083625, 63.63829205788585596792,
     24.57454391843250007810},
    {33.85296761115939290221, 0.90268363795106454361, 0.00959841694949515982, 60.43756773706766693977,
     29.25949560100913870997},
    {25.71984753303799919877, 0.89676534188092282829, 0.01062476272658970100, 61.07584814152941277143,
     28.07655374254172997439},
    {18.69037392254148954862, 0.89807977297998664579, 0.00914711061128260436, 49.71496317190430147548,
     22.35546477436284718010}};

#endif
