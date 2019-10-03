#include <iostream>
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "THn.h"
#include "TStyle.h"

TCanvas *acc_corrected(std::string data_root, std::string sim_root) {
  gStyle->SetPalette(kDarkRainBow);
  TFile *root_data = new TFile(data_root.c_str());
  TFile *root_sim = new TFile(sim_root.c_str());

  THnD *ndHist = (THnD *)root_data->Get("ndhist");
  THnD *ndHist_rec = (THnD *)root_sim->Get("ndhist");
  ndHist_rec->Scale(1 / ndHist_rec->GetEntries());
  THnD *ndHist_thrown = (THnD *)root_sim->Get("ndhist_mc");
  ndHist_thrown->Scale(1 / ndHist_thrown->GetEntries());

  for (int i = 0; i < 20; i++) ndHist_thrown->Projection(0)->SetBinContent(i, 1.0);

  THnD *ndHist_acc = (THnD *)ndHist_rec->Clone();
  ndHist_acc->Divide(ndHist_thrown);

  TCanvas *canvas = new TCanvas("c1", "c1", 1600, 900);
  canvas->Divide(2, 2);

  canvas->cd(1);
  auto W = (TH1D *)ndHist->Projection(0);
  W->Draw();

  canvas->cd(2);
  auto acc = (TH1D *)ndHist_rec->Projection(0);
  auto thrown_W = (TH1D *)ndHist_thrown->Projection(0);
  acc->Divide(thrown_W);
  acc->Draw("PMC  E1");

  canvas->cd(3);
  auto W_real = (TH1D *)ndHist->Projection(0);
  W_real->SetTitle("Acceptance Corrected");
  W_real->Multiply(acc);
  W_real->Draw("X0 E1");

  canvas->cd(4);
  auto W_scale = (TH1D *)ndHist->Projection(0);
  W_scale->Scale(1 / W_scale->GetMaximum());
  W_real->Scale(1 / W_real->GetMaximum());

  W_scale->Draw("PLC PMC PFC E1");
  W_real->Draw("SAME PLC PMC PFC E1");

  return canvas;
}

#if not defined(__CLING__)
int main(int argc, char const *argv[]) {
  if (argc < 3) {
    std::cerr << "Not enough arguments" << std::endl;
    std::cerr << "To Use:\tmc_vs_data_plots data.root mc.root" << std::endl;
    exit(1);
  }

  auto can = acc_corrected(argv[1], argv[2]);
  can->SaveAs("acc_corrected.pdf");

  return 0;
}
#endif