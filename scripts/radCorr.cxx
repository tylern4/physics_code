#include <fstream>
#include <iostream>
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TStyle.h"

static const float E1D_E0 = 4.81726;
static const float MASS_P = 0.93827203;
static const float MASS_N = 0.93956556;
static const float MASS_E = 0.000511;
static const short MAX_PARTS = 10;

// Calcuating Q^2
//	Gotten from t channel
// -q^mu^2 = -(e^mu - e^mu')^2 = Q^2
double Q2_calc(const TLorentzVector &e_mu, const TLorentzVector &e_mu_prime) {
  TLorentzVector q_mu = (e_mu - e_mu_prime);
  return -q_mu.Mag2();
}
//	Calcualting W
//	Gotten from s channel [(gamma + P)^2 == s == w^2]
//	Sqrt√[M_p^2 - Q^2 + 2 M_p gamma]
double W_calc(const TLorentzVector &e_mu, const TLorentzVector &e_mu_prime) {
  TLorentzVector q_mu = (e_mu - e_mu_prime);
  TLorentzVector p_mu(0, 0, 0, MASS_P);
  return (p_mu + q_mu).Mag();
}

int _nprt;
int _pidpart[MAX_PARTS];    //[nprt]
float _xpart[MAX_PARTS];    //[nprt]
float _ypart[MAX_PARTS];    //[nprt]
float _zpart[MAX_PARTS];    //[nprt]
float _epart[MAX_PARTS];    //[nprt]
float _pxpart[MAX_PARTS];   //[nprt]
float _pypart[MAX_PARTS];   //[nprt]
float _pzpart[MAX_PARTS];   //[nprt]
float _qpart[MAX_PARTS];    //[nprt]
int _flagspart[MAX_PARTS];  //[nprt]

void init(std::shared_ptr<TChain> _tree) {
  _tree->SetBranchStatus("*", 0);
  _tree->SetBranchAddress("epart", _epart);
  _tree->SetBranchAddress("pxpart", _pxpart);
  _tree->SetBranchAddress("pypart", _pypart);
  _tree->SetBranchAddress("pzpart", _pzpart);
}

TLorentzVector e_mu(0.0, 0.0, sqrtf(E1D_E0 *E1D_E0 - MASS_E * MASS_E), E1D_E0);
//////////////// W , Q2, Theta_star_pip, Phi_star_pip
int radCorr(const std::string &norad_root, const std::string &rad_root) {
  std::ofstream myfile;
  myfile.open("result.csv", std::ios::out | std::ios::app);
  if (!myfile.is_open()) {
    return 1;
  }

  auto noradChain = std::make_shared<TChain>("h10");
  noradChain->Add(norad_root.c_str());
  init(noradChain);

  auto wVsQ2_norad = new TH2D("wvsq2_norad", "wvsq2_norad", 500, 1.1, 2.0, 500, 1.0, 3.5);
  auto wVsQ2_rad = new TH2D("wvsq2_rad", "wvsq2_rad", 500, 1.1, 2.0, 500, 1.0, 3.5);

  auto thetaPhi_norad = new TH2D("thetaPhi_norad", "thetaPhi_norad", 500, 1.1, 2.0, 500, 1.0, 3.5);
  auto thetaPhi_rad = new TH2D("thetaPhi_rad", "thetaPhi_rad", 500, 1.1, 2.0, 500, 1.0, 3.5);

  size_t numNoRad = noradChain->GetEntries();
  for (size_t part = 0; part < numNoRad; part++) {
    if (part % 10000 == 0) std::cout << part << "\r\r" << std::flush;
    noradChain->GetEntry(part);
    TLorentzVector e_mu_prime(_pxpart[0], _pypart[0], _pzpart[0], _epart[0]);
    auto W = W_calc(e_mu, e_mu_prime);
    auto Q2 = Q2_calc(e_mu, e_mu_prime);
    wVsQ2_norad->Fill(W, Q2);
    if (std::isnan(W) || std::isnan(Q2)) continue;
    auto theta = 0;
    auto phi = 0;
    thetaPhi_norad->Fill(theta, phi);
    myfile << "norad," << W << "," << Q2 << "," << theta << "," << phi << "\n";
  }

  auto radChain = std::make_shared<TChain>("h10");
  std::cout << rad_root << std::endl;
  radChain->Add(rad_root.c_str());
  init(radChain);
  size_t numRad = radChain->GetEntries();

  if (numRad < numNoRad) {
    std::cout << numNoRad << " " << numRad << std::endl;
    return 1;
  }

  for (size_t part = 0; part < numNoRad; part++) {
    if (part % 10000 == 0) std::cout << part << "\r\r" << std::flush;
    radChain->GetEntry(part);
    TLorentzVector e_mu_prime(_pxpart[0], _pypart[0], _pzpart[0], _epart[0]);
    auto W = W_calc(e_mu, e_mu_prime);
    auto Q2 = Q2_calc(e_mu, e_mu_prime);
    wVsQ2_rad->Fill(W, Q2);
    if (std::isnan(W) || std::isnan(Q2)) continue;
    auto theta = 0;
    auto phi = 0;
    thetaPhi_rad->Fill(theta, phi);
    myfile << "rad," << W << "," << Q2 << "," << theta << "," << phi << "\n";
  }

  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  wVsQ2_norad->Draw("COLZ");
  c1->Show();

  TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
  wVsQ2_rad->Draw("COLZ");
  c2->Show();

  myfile.close();

  return 0;
}

#if not defined(__CLING__)
int main(int argc, char const *argv[]) {
  if (argc < 2) {
    std::cerr << "Not enough arguments" << std::endl;
    std::cerr << "To Use:\t" << argv[0] << " norad.root rad.root" << std::endl;
    exit(1);
  }

  return radCorr(argv[1], argv[2]);
}
#endif