#include <iostream>
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "THn.h"
#include "THnSparse.h"
#include "TStyle.h"

#define DIMENSIONS 4
//////////////// W , Q2, Theta_star_pip, Phi_star_pip
int nbins[DIMENSIONS] = {5, 5, 5, 25};
double xmin[DIMENSIONS] = {1.0, 0.8, 0, 0};
double xmax[DIMENSIONS] = {2.0, 2.0, TMath::Pi(), 2.0 * TMath::Pi()};

TCanvas *kinematics(std::string data_root) {
  TFile *root_data = new TFile(data_root.c_str());
  TCanvas *can = new TCanvas();
  can->cd();

  THnSparse *ndHist = (THnSparse *)root_data->Get("ndhist");
  // ndHist->Projection(2)->Draw();

  std::cout << "w,q2,theta,phi,y,yerr" << std::endl;
  for (int w = 0; w < nbins[0]; w++) {
    for (int q2 = 0; q2 < nbins[1]; q2++) {
      for (int theta = 0; theta < nbins[2]; theta++) {
        for (int phi = 0; phi < nbins[3]; phi++) {
          int bin[DIMENSIONS] = {w, q2, theta, phi};
          if (ndHist->GetBinContent(bin) == 0) continue;
          std::cout << w * (1.0f / nbins[0] * 1.0) + 1.0 << ",";
          std::cout << q2 * ((2.0f - 0.8f) / nbins[1]) + 0.8 << ",";
          std::cout << theta * (TMath::Pi() / nbins[2]) << ",";
          std::cout << phi * (TMath::TwoPi() / nbins[3]) << ",";
          std::cout << ndHist->GetBinContent(bin) << ",";
          std::cout << ndHist->GetBinError(ndHist->GetBin(bin)) << std::endl;
        }
      }
    }
  }

  ndHist->Projection(3)->Draw("COLZ");

  return can;
}

#if not defined(__CLING__)
int main(int argc, char const *argv[]) {
  if (argc < 2) {
    std::cerr << "Not enough arguments" << std::endl;
    std::cerr << "To Use:\t" << argv[0] << " data.root" << std::endl;
    exit(1);
  }

  auto can = kinematics(argv[1]);
  can->SaveAs("kinematics.pdf");

  return 0;
}
#endif