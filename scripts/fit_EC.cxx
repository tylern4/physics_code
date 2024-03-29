#include <iostream>
#include "Math/Minimizer.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TLine.h"
#include "TROOT.h"
#include "TStyle.h"

#define THETA_BINS 100
#define PHI_BINS 100

#define NUM_SECTORS 6
#define FID_SLICES 40

const double min_phi[6] = {60, 0, -60, -120, -180, 120};
const double max_phi[6] = {120, 60, 0, -60, -120, 180};
TFile *out = new TFile("out.root", "RECREATE");

struct FitsData {
  double right_edge;
  double left_edge;
  TF1 *fit;
};

double genNormal(double *x, double *par) {
  double frac = (par[1] / (2 * par[0] * TMath::Gamma(1 / par[1])));
  double expo = TMath::Power(TMath::Abs(x[0] - par[2]) / par[0], par[1]);

  double func = par[3] * frac * TMath::Exp(-expo);

  return func;
}

double gaus(double *x, double *par) {
  double frac = (par[1] / (2 * par[0] * TMath::Gamma(1 / par[1])));
  double expo = TMath::Power(TMath::Abs(x[0] - par[2]) / par[0], par[1]);

  double func = par[3] * frac * TMath::Exp(-expo);

  return func;
}

FitsData FitGenNormal(TH1D *hist, int slice, int sec_i) {
  double min_value = hist->GetXaxis()->GetXmin();
  double max_value = hist->GetXaxis()->GetXmax();

  TCanvas *can =
      new TCanvas(Form("can_sec%d_slice%d", sec_i, slice), Form("can_sec%d_slice%d", sec_i, slice), 800, 600);
  out->cd();
  can->cd();

  hist->Rebin(2);
  TF1 *fitFunc = new TF1("genNormal", genNormal, min_value, max_value, 4);
  double min, max, val, min_m, max_m, min_fwhm, max_fwhm;

  fitFunc->SetParLimits(1, 2.0, 200.0);
  fitFunc->SetParameter(0, 15.0);
  fitFunc->SetParameter(1, 10.0);
  fitFunc->SetParameter(2, (min_value + max_value) / 2.0);
  fitFunc->SetParameter(3, 3000.0);
  fitFunc->SetParNames("alpha", "beta", "mu", "weight");

  for (int i = 0; i < 10; i++) hist->Fit("genNormal", "QMR0+", "", min_value, max_value);

  hist->Fit("genNormal", "QMR+", "", min_value, max_value);
  double fwhm_max = fitFunc->GetMaximum() / 2;

  min_fwhm = fitFunc->GetX(fitFunc->GetMaximum() / 2.0, min_value, max_value / 2.0, 1.E-10, 100000, false);
  max_fwhm = fitFunc->GetX(fitFunc->GetMaximum() / 2.0, max_value / 2.0, max_value, 1.E-10, 100000, false);

  for (double m = min_value; m < max_value; m = m + 0.0001) {
    val = fitFunc->Derivative(m);
    if (max < val) {
      max = val;
      max_m = m;
    }
    if (val < min) {
      min = val;
      min_m = m;
    }
  }

  float ymax = hist->GetMaximum();

  TLine *lineR = new TLine(max_m, 0, max_m, ymax);
  lineR->SetLineColor(kBlue);
  lineR->SetLineWidth(3);

  TLine *lineL = new TLine(min_m, 0, min_m, ymax);
  lineL->SetLineColor(kBlue);
  lineL->SetLineWidth(3);

  TLine *lineR_fwhm = new TLine(max_fwhm, 0, max_fwhm, ymax);
  lineR_fwhm->SetLineColor(kGreen);
  lineR_fwhm->SetLineWidth(3);

  TLine *lineL_fwhm = new TLine(min_fwhm, 0, min_fwhm, ymax);
  lineL_fwhm->SetLineColor(kGreen);
  lineL_fwhm->SetLineWidth(3);

  hist->Draw();
  fitFunc->Draw("same");
  lineL->Draw("same");
  lineR->Draw("same");
  lineL_fwhm->Draw("same");
  lineR_fwhm->Draw("same");
  can->Write();

  FitsData output;
  output.left_edge = min_m;
  output.right_edge = max_m;
  output.fit = fitFunc;
  return output;
}

double poly(double *x, double *par) {
  double func = 0.0;
  func += -par[2];
  func += -(TMath::Log(par[3] * (1 - (x[0] / par[0])))) / par[1];
  return func;
}

double _fid(double *x, double *par) {
  double func = 0.0;
  func += par[0] * (1 + par[3] * TMath::Exp(-par[1] * (x[0] + par[2])));

  return func;
}

FitsData FitTheFit(TGraph *hist, int sec, bool lr) {
  double min_value = min_phi[sec];
  double max_value = max_phi[sec];
  double center = (min_phi[sec] + max_phi[sec]) / 2.0;

  if (lr) {
    min_value = min_phi[sec];
    max_value = (min_value + max_value) / 2.0;
  } else {
    min_value = (min_value + max_value) / 2.0;
    max_value = max_phi[sec];
  }
  min_value = 30;
  max_value = 60;

  TF1 *fitFunc = new TF1("fid", _fid, min_value, max_value, 4);
  if (lr) {
    fitFunc->SetLineColor(kGreen);
    fitFunc->FixParameter(3, -1);
  } else {
    fitFunc->FixParameter(3, 1);
  }
  double min, max, val, min_m, max_m;

  fitFunc->SetParNames("a", "b", "c", "d", "e");

  for (int i = 0; i < 10; i++) hist->Fit("fid", "QMR0+", "", min_value, max_value);

  hist->Fit("fid", "QMR+", "", min_value, max_value);

  FitsData output;
  output.left_edge = min_m;
  output.right_edge = max_m;
  output.fit = fitFunc;
  return output;
}

int fit_EC(std::string file = "e1d_all.root") {
  TFile *f = new TFile(file.c_str());
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(1111);
  gStyle->SetNumberContours(90);
  gStyle->SetPadGridX(kTRUE);
  gStyle->SetPadGridY(kTRUE);
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");

  TH2D *elec_fid[NUM_SECTORS];
  for (int i = 0; i < NUM_SECTORS; i++) {
    elec_fid[i] = (TH2D *)f->Get(Form("Fid_cuts/electron_fid_sec%d", i + 1));
  }
  int bin_nums = 500;
  double slice_width = ((double)bin_nums / (double)FID_SLICES);
  double y_width = (60.0 / (double)bin_nums);
  TH1D *electron_fid_sec_slice[NUM_SECTORS][FID_SLICES];
  TH1D *slices[4][NUM_SECTORS];
  TGraph *fid_graph[NUM_SECTORS];

  for (int sec_i = 0; sec_i < NUM_SECTORS; sec_i++) {
    std::vector<double> x;
    std::vector<double> y;

    for (int slice = 14; slice < FID_SLICES; slice++) {
      electron_fid_sec_slice[sec_i][slice] =
          (TH1D *)elec_fid[sec_i]->ProjectionX(Form("electron_fid_sec_%d_%d", sec_i + 1, slice + 1),
                                               slice_width * slice, slice_width * slice + (slice_width - 1));
      if (electron_fid_sec_slice[sec_i][slice]->GetEntries() > 100) {
        FitsData genNormal_fit = FitGenNormal(electron_fid_sec_slice[sec_i][slice], slice, sec_i);

        x.push_back(y_width * slice_width * slice);
        y.push_back(genNormal_fit.left_edge);
        x.push_back(y_width * slice_width * slice);
        y.push_back(genNormal_fit.right_edge);
      }
    }
    fid_graph[sec_i] = new TGraph(y.size(), &x[0], &y[0]);
  }

  TCanvas *fid_can = new TCanvas("fid_can", "fid_can", 1600, 900);
  fid_can->DivideSquare(6);

  for (int i = 0; i < NUM_SECTORS; i++) {
    fid_can->cd(i + 1);
    elec_fid[i]->Draw("colz");
    fid_graph[i]->Draw("same*");
    FitTheFit(fid_graph[i], i, true).fit->Draw("same");
    FitTheFit(fid_graph[i], i, false).fit->Draw("same");
  }
  fid_can->Write();
  out->Write();
  // gROOT->ProcessLine(".q;");
  return 0;
}

#if not defined(__CLING__)
int main(int argc, char const *argv[]) {
  if (argc < 2) {
    std::cerr << "Not enough arguments" << std::endl;
    std::cerr << "To Use:\t" << argv[0] << " v2_all.root" << std::endl;
    exit(1);
  }

  fit_EC(argv[1]);

  return 0;
}
#endif