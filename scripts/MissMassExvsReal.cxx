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

TCanvas *MissMassExvsReal(const std::string &data_root, const std::string &mc_root) {
  TFile *root_data = new TFile(data_root.c_str());
  TFile *root_mc = new TFile(mc_root.c_str());
  TCanvas *can = new TCanvas();
  can->Divide(3, 2);
  for (int i = 0; i < 6; i++) {
    can->cd(i + 1);
    TH1D *data = (TH1D *)root_data->Get(Form("Missing_Mass/Missing_Mass_small_%d", i));
    TH1D *mc = (TH1D *)root_mc->Get(Form("Missing_Mass/Missing_Mass_small_%d", i));

    data->Scale(1.0 / data->GetMaximum());
    mc->Scale(1.0 / mc->GetMaximum());
    mc->SetLineColor(kRed);
    mc->SetMarkerColor(kRed);

    data->Draw();
    mc->Draw("SAME");
  }

  return can;
}

#if not defined(__CLING__)
int main(int argc, char const *argv[]) {
  if (argc < 2) {
    std::cerr << "Not enough arguments" << std::endl;
    std::cerr << "To Use:\t" << argv[0] << " data.root mc.root" << std::endl;
    exit(1);
  }

  auto can = MissMassExvsReal(argv[1], argv[2]);
  can->SaveAs("MissMassExvsReal.pdf");

  return 0;
}
#endif