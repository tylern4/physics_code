#include <iostream>
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "THn.h"
#include "THnSparse.h"
#include "TStyle.h"

//////////////// W , Q2, Theta_star_pip, Phi_star_pip

TCanvas *kinematics(std::string data_root) {
  TFile *root_data = new TFile(data_root.c_str());
  TCanvas *can = new TCanvas();
  can->cd();

  THnSparse *ndHist = (THnSparse *)root_data->Get("ndhist");
  const int DIMENSIONS = ndHist->GetNdimensions();
  int nbins[DIMENSIONS];
  double xmin[DIMENSIONS];
  double xmax[DIMENSIONS];
  int bin[DIMENSIONS];
  for (int i = 0; i < DIMENSIONS; i++) nbins[i] = ndHist->GetAxis(i)->GetNbins();
  for (int i = 0; i < DIMENSIONS; i++) xmin[i] = ndHist->GetAxis(i)->GetBinLowEdge(1);
  for (int i = 0; i < DIMENSIONS; i++) xmax[i] = ndHist->GetAxis(i)->GetBinUpEdge(ndHist->GetAxis(i)->GetNbins());

  ndHist->Projection(2)->Draw();

  std::cout << "q2,w,theta,phi,y,yerr" << std::endl;
  for (int w = 0; w < nbins[0]; w++) {
    for (int q2 = 0; q2 < nbins[1]; q2++) {
      for (int theta = 0; theta < nbins[2]; theta++) {
        for (int phi = 0; phi < nbins[3]; phi++) {
          bin[0] = w;
          bin[1] = q2;
          bin[2] = theta;
          bin[3] = phi;
          if (ndHist->GetBinContent(bin) == 0) continue;
          std::cout << q2 * ((xmax[1] - xmin[1]) / nbins[1] * 1.0) + xmin[1] << ",";
          std::cout << w * ((xmax[0] - xmin[0]) / nbins[0] * 1.0) + xmin[0] << ",";
          std::cout << theta * ((xmax[2] - xmin[2]) / nbins[2] * 1.0) + xmin[2] << ",";
          std::cout << phi * ((xmax[3] - xmin[3]) / nbins[3] * 1.0) + xmin[3] << ",";
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