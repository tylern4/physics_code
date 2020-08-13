#include <iostream>
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "THn.h"
#include "THnSparse.h"
#include "TStyle.h"

TCanvas *MomCorrection(std::string with, std::string without) {
  TFile *with_data = new TFile(with.c_str());
  TFile *without_data = new TFile(without.c_str());

  auto can = new TCanvas("mom_corrections", "mom_corrections", 1600, 900);
  can->Divide(3, 2);
  for (int sec = 1; sec <= 6; sec++) {
    can->cd(sec);
    TH1 *temp1 = (TH1 *)with_data->Get(Form("W_vs_Q2_sec/W_sec_%d", sec));
    temp1->Scale(1 / temp1->GetMaximum());
    temp1->SetLineColor(kRed);
    temp1->Draw();
    TH1 *temp2 = (TH1 *)without_data->Get(Form("W_vs_Q2_sec/W_sec_%d", sec));
    temp2->Scale(1 / temp2->GetMaximum());
    temp2->Draw("same");
  }

  auto can1 = new TCanvas("MM_mom_corrections", "MM_mom_corrections", 1600, 900);
  can1->Divide(3, 2);
  for (int sec = 1; sec <= 6; sec++) {
    can1->cd(sec);
    TH1 *temp1 = (TH1 *)with_data->Get(Form("Missing_Mass/Missing_Mass_small_%d", sec - 1));
    temp1->Scale(1 / temp1->GetMaximum());
    temp1->SetLineColor(kRed);
    temp1->Draw();
    TH1 *temp2 = (TH1 *)without_data->Get(Form("Missing_Mass/Missing_Mass_small_%d", sec - 1));
    temp2->Scale(1 / temp2->GetMaximum());
    temp2->Draw("same");
  }

  auto can2 = new TCanvas("channel_mom_corrections", "channel_mom_corrections", 1600, 900);
  can2->Divide(3, 2);
  for (int sec = 1; sec <= 6; sec++) {
    can2->cd(sec);
    TH1 *temp1 = (TH1 *)with_data->Get(Form("W_vs_Q2_sec/W_channel_sec_%d", sec));
    temp1->Scale(1 / temp1->GetMaximum());
    temp1->SetLineColor(kRed);
    temp1->GetXaxis()->SetRangeUser(1.2, 2.0);
    temp1->Draw();
    TH1 *temp2 = (TH1 *)without_data->Get(Form("W_vs_Q2_sec/W_channel_sec_%d", sec));
    temp2->Scale(1 / temp2->GetMaximum());
    temp2->GetXaxis()->SetRangeUser(1.2, 2.0);
    temp2->Draw("same");
  }

  return can;
}

#if not defined(__CLING__)
int main(int argc, char const *argv[]) {
  if (argc < 3) {
    std::cerr << "Not enough arguments" << std::endl;
    std::cerr << "To Use:\t" << argv[0] << " with_momCorr.root no_momCorr.root" << std::endl;
    exit(1);
  }

  auto can = MomCorrection(argv[1], argv[2]);
  can->SaveAs("MomCorrection.pdf");

  return 0;
}
#endif