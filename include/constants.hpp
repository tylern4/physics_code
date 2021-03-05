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

#define PRINT_TIMEING(start_time, name)                                                                             \
  std::cout                                                                                                         \
      << name                                                                                                       \
      << static_cast<std::chrono::duration<double>>(std::chrono::high_resolution_clock::now() - start_time).count() \
      << " sec" << std::endl

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<float>> LorentzVector;
typedef ROOT::Math::Rotation3D LorentzRotation;
typedef ROOT::Math::XYZVectorD Vector3;

static const short NUM_THREADS = (getenv("NUM_THREADS") != NULL) ? atoi(getenv("NUM_THREADS")) : 8;
// static const float BEAM_E = (getenv("BEAM_E") != NULL) ? atof(getenv("BEAM_E")) : 4.81726;

static const short NUM_SECTORS = 6;
static const short MAX_PARTS = 100;
static const short N_SIGMA = 3;
static const float PI = TMath::Pi();
static const float INV_SQRT_2PI = TMath::Sqrt(2 * TMath::Pi());
static const float D2R = PI / 180.0;
static const float DEG2RAD = PI / 180.0;
static const float RAD2DEG = 180.0 / PI;
static const int POSITIVE = 1;
static const int NEGATIVE = -1;
// static const float E1D_E0 = 4.802;   // GeV
// Run Info 22848 - 23500
// http://clasweb.jlab.org/clasonline/servlet/runinfo?action=detail&run=22880
static const float E1D_E0 = 4.81726;  // GeV
static const float E1F_E0 = 5.479;    // GeV
static const float E16_E0 = 5.76959;

static const float SOL = 29.9792458;
// misc. constants
static const float FSC = 0.00729735253;
static const float NA = 6.02214129E23;               // Avigadro's number
static const float QE = 1.60217646E-19;              // Charge or electron
static const float FS_ALPHA = 0.007297352570866302;  // Fine structure alpha

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
static const float MASS_P = 0.93827203;
static const float MASS_N = 0.93956556;
static const float MASS_E = 0.000511;
static const float MASS_PIP = 0.13957018;
static const float MASS_PIM = 0.13957018;
static const float MASS_PI0 = 0.1349766;
static const float MASS_KP = 0.493677;
static const float MASS_KM = 0.493677;
static const float MASS_G = 0.0;
static const float MASS_OMEGA = 0.78265;

// Got from ye/arjun
// TODO::Get right EC threshhold
// p_min(in MeV) = 214 + 2.47 × EC_threshold(in mV)
static const float MIN_P_CUT = 0.005f;

// const float min_phi[6] = {0, 60, 120, -180, -120, -60};
// const float max_phi[6] = {60, 120, 180, -120, -60, 0};
enum fid_part { hadron, proton, pip };

static std::unordered_map<int, std::pair<float, float>> phi_map = {
    {1, std::make_pair(60, 120)},   {2, std::make_pair(0, 60)},      {3, std::make_pair(-60, 0)},
    {4, std::make_pair(-120, -60)}, {5, std::make_pair(-180, -120)}, {6, std::make_pair(120, 180)}};

static std::unordered_map<int, float> phi_center = {{1, 90}, {2, 30}, {3, -30}, {4, -90}, {5, -150}, {6, 150}};

const float min_phi[6] = {60, 0, -60, -120, -180, 120};
const float max_phi[6] = {120, 60, 0, -60, -120, 180};

static std::unordered_map<int, float> mass_map = {
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

static const std::vector<double> dt_P_top = {-2.274399517322955, 2.302723100001534, -0.318125383189543,
                                             2.4275805364424268, 0.7291582349871679};
static const std::vector<double> dt_P_bot = {3.1816135800162315, 1.1925931756674077, 0.4959408417527087,
                                             -3.5510535607301335, 1.4340601178263646};

// log_sqrt_pol1
// static const std::vector<double> dt_pip_top = {-0.5680719469021316, 4.065155482173403,  -2.285775227566277,
//                                                0.09291303310648173, 1.6698856065697703, 0.6181398119277808};

// log_pol2
static const std::vector<double> dt_pip_top = {-0.8735623149163128, 0.13773880009866746, -0.26618769381591156,
                                               1.4340820263080105, -1.9504462977780122};
// log_sqrt_pol1
// static const std::vector<double> dt_pip_bottom = {0.7021097407492168, 1.1754282714198145,  0.851129051933132,
//                                                   0.8804059391935983, -1.9632091423474742, 0.22761907815192178};

// log_pol2
static const std::vector<double> dt_pip_bottom = {0.8668333301141518, 0.8935201118434752, 0.2310903343362912,
                                                  -1.3830410597583622, 0.26273043482794134};

enum mm2_names { A, mu, sigma, lambda1, lambda2 };
static const std::vector<std::vector<float>> mm2_fit_values = {
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

// static const std::array< mom_corr_electron[6][4][3];
// static const std::array<std::array<double, 3>, 3> arr = {{{5, 8, 2}, {8, 3, 1}, {5, 3, 9}}};

static const double mom_corr_electron[6][5][3] = {
    {{-1.2001951889148072e-07, 4.36084590672119e-06, -3.954321133147653e-05},
     {4.440732524207208e-07, -1.5545573051829904e-05, 0.00013544408606948316},
     {9.629921223279887e-07, -3.580761848014096e-05, 0.0003312867905874558},
     {-5.65507875183858e-06, 0.00020216769804674918, -0.0018148689572218206},
     {1.5769436445717492e-05, -0.0005625490036981456, 0.00426151889051263}},
    {{-4.815126880139836e-07, 1.7000832091845205e-05, -0.0001495256121353645},
     {9.457558735006025e-07, -3.356404644193405e-05, 0.0002965694884907279},
     {5.7499896735362e-06, -0.00020069724889433608, 0.001743174929218162},
     {-6.710655239052169e-06, 0.0002497990095547633, -0.002295624539162647},
     {3.816177587969378e-05, -0.0014144900869375674, 0.012669636262929967}},
    {{-3.3778845138411547e-07, 1.1986123146786749e-05, -0.00010596768947051746},
     {-1.094199590711473e-06, 3.8919602174231565e-05, -0.00034489508368228746},
     {4.344329871751627e-06, -0.0001542279916664165, 0.0013625502584535731},
     {2.353425938104113e-06, -9.3354891510259e-05, 0.001012781455622935},
     {2.530790018490971e-06, -0.00016298994883765438, 0.0012304011128389822}},
    {{-1.402626910706675e-07, 5.072163392549302e-06, -4.575808992782373e-05},
     {-1.2077650843461844e-06, 4.2841458004755126e-05, -0.00037855427564291325},
     {2.4512750839939776e-06, -8.988061343143945e-05, 0.0008210157727153251},
     {9.111330595592564e-06, -0.0003224828506224554, 0.0029204401379112756},
     {4.212026024777228e-05, -0.0014358842628037648, 0.01139813725587655}},
    {{6.613315299555605e-07, -2.342031410026976e-05, 0.0002066124279934897},
     {-9.695867296910675e-07, 3.42539589290092e-05, -0.0003013865348151469},
     {-7.590024168948459e-06, 0.00026897192445843505, -0.0023733343964052827},
     {3.0029302689265534e-06, -0.00010155991725819358, 0.0008611221253852591},
     {2.6965260005103033e-05, -0.0008690561301267078, 0.005753104854809957}},
    {{2.3204189129383225e-08, -9.579986465422747e-07, 9.75317191426648e-06},
     {-6.164779266373761e-07, 2.1681375099722985e-05, -0.0001899482281826151},
     {2.0455924843919723e-06, -6.779754801886723e-05, 0.0005548703649520899},
     {1.1294969699847741e-06, -3.505495431321111e-05, 0.00032422156947394556},
     {8.067053769649876e-06, -0.00030569606835432935, 0.00169506129639814}},

};

#endif
