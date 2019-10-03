#include <iostream>
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "THn.h"
#include "TStyle.h"

TCanvas *kinematics(std::string data_root, std::string sim_root) {
  // gStyle->SetPalette(kDarkRainBow);
  TFile *root_data = new TFile(data_root.c_str());

  TCanvas *canvas = new TCanvas("c1", "c1", 1600, 900);
  canvas->cd();

  THnD *ndHist = (THnD *)root_data->Get("ndhist");
  ndHist->Projection(1, 0)->Draw("COLZ");

  return canvas;
}

#if not defined(__CLING__)
int main(int argc, char const *argv[]) {
  if (argc < 3) {
    std::cerr << "Not enough arguments" << std::endl;
    std::cerr << "To Use:\tmc_vs_data_plots data.root mc.root" << std::endl;
    exit(1);
  }

  auto can = kinematics(argv[1], argv[2]);
  can->SaveAs("kinematics.pdf");

  return 0;
}
#endif