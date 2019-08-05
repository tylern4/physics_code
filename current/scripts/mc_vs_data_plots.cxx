#include <iostream>
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"

TCanvas *mc_vs_data_plots(std::string data_root, std::string mc_root) {
  gStyle->SetPalette(kDarkRainBow);
  TFile *root_data = new TFile(data_root.c_str());
  TFile *root_mc = new TFile(mc_root.c_str());
  TCanvas *canvas = new TCanvas("c1", "c1", 1600, 900);
  TH1D *W_data = (TH1D *)root_data->Get("W vs Q2/W_channel");
  TH1D *W_mc_data = (TH1D *)root_mc->Get("W vs Q2/W_channel");
  TH1D *W_mc_true = (TH1D *)root_mc->Get("W vs Q2 MC/W_MC");

  W_data->GetXaxis()->SetRangeUser(1.1, 1.8);

  Double_t norm = 1;
  W_data->Scale(norm / W_data->GetMaximum());
  W_mc_data->Scale(norm / W_mc_data->GetMaximum());
  W_mc_true->Scale(norm / W_mc_true->GetMaximum());

  W_data->SetMarkerStyle(kFullCircle);
  W_mc_true->SetMarkerStyle(kFullCircle);
  W_mc_data->SetMarkerStyle(kFullCircle);

  W_data->Draw("PLC PMC");
  W_mc_data->Draw("SAME PLC PMC");
  W_mc_true->Draw("SAME PLC PMC");

  return canvas;
}

#if !defined(__CLING__)
int main(int argc, char const *argv[]) {
  if (argc < 3) {
    std::cerr << "Not enough arguments" << std::endl;
    std::cerr << "To Use:\tmc_vs_data_plots data.root mc.root" << std::endl;
    exit(1);
  }

  auto can = mc_vs_data_plots(argv[1], argv[2]);
  can->SaveAs("mc_vs_data_plots.pdf");
  return 0;
}
#endif