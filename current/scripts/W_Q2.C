#include <iostream>
#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"

int W_Q2() {
  TFile *f_wFid = new TFile("v2_all.root");
  TFile *f_noFid = new TFile("v2_all_noFid.root");

  TCanvas *wq2 = new TCanvas("wq2", "wq2", 1600, 900);
  wq2->Divide(0, 3);

  TH2D *wq2_noFid = (TH2D *)f_noFid->Get("W vs Q2/WvsQ2_channel");
  TH2D *wq2_wFid = (TH2D *)f_wFid->Get("W vs Q2/WvsQ2_channel");
  TH2D *wq2_sub = (TH2D *)wq2_noFid->Clone();
  wq2_sub->Add(wq2_wFid, -1);

  wq2->cd(1);
  wq2_noFid->Draw("same");
  wq2->cd(2);
  wq2_wFid->Draw("same");
  wq2->cd(3);
  wq2_sub->Draw("same");

  TCanvas *w_can = new TCanvas("w_can", "w_can", 1600, 900);
  w_can->Divide(0, 3);

  TH1D *w_noFid = (TH1D *)f_noFid->Get("W vs Q2/W_channel");
  TH1D *w_wFid = (TH1D *)f_wFid->Get("W vs Q2/W_channel");
  TH1D *w_sub = (TH1D *)w_noFid->Clone();
  w_sub->Add(w_wFid, -1);

  w_can->cd(1);
  w_noFid->Draw("same");
  w_can->cd(2);
  w_wFid->Draw();
  w_can->cd(3);
  w_noFid->SetLineColor(kBlack);
  w_noFid->Draw("same");
  w_wFid->SetLineColor(kBlue);
  w_wFid->Draw("same");
  w_sub->SetLineColor(kRed);
  w_sub->Draw("same");

  TCanvas *q2_can = new TCanvas("q2_can", "q2_can", 1600, 900);
  q2_can->Divide(0, 3);

  TH1D *q2_noFid = (TH1D *)f_noFid->Get("W vs Q2/Q2_channel");
  TH1D *q2_wFid = (TH1D *)f_wFid->Get("W vs Q2/Q2_channel");
  TH1D *q2_sub = (TH1D *)q2_noFid->Clone();
  q2_sub->Add(q2_wFid, -1);

  q2_can->cd(1);
  q2_noFid->Draw("same");
  q2_can->cd(2);
  q2_wFid->Draw();
  q2_can->cd(3);
  q2_noFid->SetLineColor(kBlack);
  q2_noFid->Draw("same");
  q2_wFid->SetLineColor(kBlue);
  q2_wFid->Draw("same");
  q2_sub->SetLineColor(kRed);
  q2_sub->Draw("same");

  TCanvas *c1 = new TCanvas("c1", "c1", 1600, 900);

  c1->Divide(0, 3);

  TH2D *fid_noFid = (TH2D *)f_noFid->Get("Fid_cuts/electron_fid");
  TH2D *fid_wFid = (TH2D *)f_wFid->Get("Fid_cuts/electron_fid");
  TH2D *fid_sub = (TH2D *)fid_noFid->Clone(); 
  fid_sub->Add(fid_wFid, -1);

  c1->cd(1);
  fid_noFid->Draw("same");
  c1->cd(2);
  fid_wFid->Draw();
  c1->cd(3);
  fid_sub->Draw();


  return 0;
}

int binned() {
  TFile *f_wFid = new TFile("v2_all.root");
  TFile *f_noFid = new TFile("v2_all_noFid.root");

  TCanvas *wq2 = new TCanvas("wq2_binned", "wq2_binned", 1600, 900);
  wq2->Divide(0, 3);

  TH2D *wq2_noFid = (TH2D *)f_noFid->Get("W_Q2_binned/WvsQ2_hist_binned");
  TH2D *wq2_wFid = (TH2D *)f_wFid->Get("W_Q2_binned/WvsQ2_hist_binned");
  TH2D *wq2_sub = (TH2D *)wq2_noFid->Clone();
  wq2_sub->Add(wq2_wFid, -1);

  wq2->cd(1);
  wq2_noFid->Draw("same");
  wq2->cd(2);
  wq2_wFid->Draw("same");
  wq2->cd(3);
  wq2_sub->Draw("same");

  return 0;
}
