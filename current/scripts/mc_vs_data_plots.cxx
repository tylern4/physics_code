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

TCanvas *plot_acceptance_corrected(std::string data_root, std::string mc_root) {
  gStyle->SetPalette(kDarkRainBow);
  TCanvas *canvas = new TCanvas("c2", "c2", 2600, 1900);
  canvas->Divide(2, 2);

  TFile *root_data = new TFile(data_root.c_str());
  TFile *root_mc = new TFile(mc_root.c_str());

  TH1D *W_data = (TH1D *)root_data->Get("W vs Q2/W_channel");
  TH1D *W_channel = (TH1D *)W_data->Clone("W_channel_data");
  TH1D *W_mc_rec = (TH1D *)root_mc->Get("W vs Q2/W_channel");
  TH1D *W_mc_thrown = (TH1D *)root_mc->Get("W vs Q2 MC/W_MC");

  Double_t norm = 1.0;
  // W_channel->Scale(norm / W_channel->GetMaximum());
  // W_data->Scale(norm / W_data->GetMaximum());
  // W_mc_rec->Scale(norm / W_mc_rec->GetMaximum());
  // W_mc_thrown->Scale(norm / W_mc_thrown->GetMaximum());

  W_channel->GetXaxis()->SetRangeUser(1.2, 1.8);
  W_data->GetXaxis()->SetRangeUser(1.2, 1.8);
  W_mc_rec->GetXaxis()->SetRangeUser(1.2, 1.8);
  W_mc_thrown->GetXaxis()->SetRangeUser(1.2, 1.8);

  canvas->cd(1);
  W_channel->Draw("SAME PLC PMC");

  TH1D *acceptance = (TH1D *)W_mc_thrown->Clone("Acceptance");
  acceptance->Divide(W_mc_rec);
  canvas->cd(2);
  acceptance->Draw("SAME PLC PMC");

  W_data->Multiply(acceptance);
  canvas->cd(3);
  W_data->Draw("SAME PLC PMC");

  canvas->cd(4);

  return canvas;
}

#if not defined(__CLING__)
int main(int argc, char const *argv[]) {
  if (argc < 3) {
    std::cerr << "Not enough arguments" << std::endl;
    std::cerr << "To Use:\tmc_vs_data_plots data.root mc.root" << std::endl;
    exit(1);
  }

  auto can = mc_vs_data_plots(argv[1], argv[2]);
  can->SaveAs("mc_vs_data_plots.pdf");

  auto can1 = plot_acceptance_corrected(argv[1], argv[2]);
  can1->SaveAs("plot_acceptance_corrected.pdf");

  return 0;
}
#endif