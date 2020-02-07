#include <iostream>
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "THn.h"
#include "TStyle.h"

TCanvas *DrawProjection(THnD *ndHist, THnD *ndHist_rec, THnD *ndHist_thrown, int projection) {
  TCanvas *canvas = new TCanvas("c1", "c1", 1600, 900);
  canvas->Divide(2, 2);
  canvas->cd(1);
  auto hist = (TH1D *)ndHist->Projection(projection);
  hist->Draw();
  canvas->cd(2);
  auto acc = (TH1D *)ndHist_rec->Projection(projection);
  auto thrown_hist = (TH1D *)ndHist_thrown->Projection(projection);
  acc->Divide(thrown_hist);
  acc->Draw("PMC  E1");
  canvas->cd(3);
  auto hist_real = (TH1D *)ndHist->Projection(projection);
  hist_real->SetTitle("Acceptance Corrected");
  hist_real->Multiply(acc);
  hist_real->Draw("X0 E1");
  canvas->cd(4);
  auto hist_scale = (TH1D *)ndHist->Projection(projection);
  hist_scale->Scale(1 / hist_scale->GetMaximum());
  hist_real->Scale(1 / hist_real->GetMaximum());
  hist_scale->Draw("PLC PMC PFC E1");
  hist_real->Draw("SAME PLC PMC PFC E1");

  return canvas;
}

TCanvas *acc_corrected(std::string data_root, std::string sim_root, int projection = 0) {
  // gStyle->SetPalette(kDarkRainBow);
  TFile *root_data = new TFile(data_root.c_str());
  TFile *root_sim = new TFile(sim_root.c_str());

  THnD *ndHist = (THnD *)root_data->Get("ndhist");
  THnD *ndHist_rec = (THnD *)root_sim->Get("ndhist");
  ndHist_rec->Scale(1 / ndHist_rec->GetEntries());
  THnD *ndHist_thrown = (THnD *)root_sim->Get("ndhist_mc");
  ndHist_thrown->Scale(1 / ndHist_thrown->GetEntries());

  return DrawProjection(ndHist, ndHist_rec, ndHist_thrown, projection);
}

TCanvas *DrawProjection_2D(THnD *ndHist, THnD *ndHist_rec, THnD *ndHist_thrown, std::pair<int, int> projection) {
  TCanvas *canvas = new TCanvas("c1", "c1", 1600, 900);
  canvas->Divide(2, 2);
  canvas->cd(1);
  auto hist = (TH2D *)ndHist->Projection(projection.first, projection.second);
  hist->Draw("COLZ");
  canvas->cd(2);
  auto acc = (TH2D *)ndHist_rec->Projection(projection.first, projection.second);
  auto thrown_hist = (TH2D *)ndHist_thrown->Projection(projection.first, projection.second);
  acc->Divide(thrown_hist);
  acc->Draw("COLZ");
  canvas->cd(3);
  auto hist_real = (TH2D *)ndHist->Projection(projection.first, projection.second);
  hist_real->SetTitle("Acceptance Corrected");
  hist_real->Multiply(acc);
  hist_real->Draw("COLZ");
  canvas->cd(4);
  auto hist_scale = (TH2D *)ndHist->Projection(projection.first, projection.second);
  hist_scale->Scale(1 / hist_scale->GetMaximum());
  hist_real->Scale(1 / hist_real->GetMaximum());
  hist_scale->Draw("COLZ");

  return canvas;
}

TCanvas *twoD_acc_corrected(std::string data_root, std::string sim_root, std::pair<int, int> projection) {
  // gStyle->SetPalette(kDarkRainBow);
  TFile *root_data = new TFile(data_root.c_str());
  TFile *root_sim = new TFile(sim_root.c_str());

  THnD *ndHist = (THnD *)root_data->Get("ndhist");
  THnD *ndHist_rec = (THnD *)root_sim->Get("ndhist");
  ndHist_rec->Scale(1 / ndHist_rec->GetEntries());
  THnD *ndHist_thrown = (THnD *)root_sim->Get("ndhist_mc");
  ndHist_thrown->Scale(1 / ndHist_thrown->GetEntries());

  ndHist->Sumw2();
  ndHist_rec->Sumw2();
  ndHist_thrown->Sumw2();

  return DrawProjection_2D(ndHist, ndHist_rec, ndHist_thrown, projection);
}

#if not defined(__CLING__)
int main(int argc, char const *argv[]) {
  if (argc < 3) {
    std::cerr << "Not enough arguments" << std::endl;
    std::cerr << "To Use:\tmc_vs_data_plots data.root mc.root" << std::endl;
    exit(1);
  }

  std::map<int, std::string> projections = {{0, "W"}, {1, "Q2"}, {2, "Theta"}, {3, "Phi"}};
  for (auto &p : projections) {
    auto can = acc_corrected(argv[1], argv[2], p.first);
    can->SaveAs(Form("acc_corrected_%s.pdf", p.second.c_str()));
    delete can;
  }

  std::map<std::pair<int, int>, std::string> projections_2d = {
      {std::make_pair(1, 0), "W vs Q^2"}, {std::make_pair(1, 2), "W vs Theta"}, {std::make_pair(2, 3), "Theta vs Phi"}};

  for (auto &p : projections_2d) {
    auto can = twoD_acc_corrected(argv[1], argv[2], p.first);
    can->SaveAs(Form("twoD_acc_corrected_%s.pdf", p.second.c_str()));
    delete can;
  }

  return 0;
}
#endif