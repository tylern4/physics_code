#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TLine.h"
#include "TROOT.h"

// Quadratic background function
Double_t background(Double_t *x, Double_t *par) {
  return par[0] + par[1] * x[0] + par[2] * x[0] * x[0];
}

// Lorenzian Peak function
Double_t lorentzianPeak(Double_t *x, Double_t *par) {
  return (0.5 * par[0] * par[1] / TMath::Pi()) /
         TMath::Max(1.e-10,
                    (x[0] - par[2]) * (x[0] - par[2]) + .25 * par[1] * par[1]);
}

// Sum of background and peak function
Double_t fitFunction(Double_t *x, Double_t *par) {
  return background(x, par) + lorentzianPeak(x, &par[3]);
}

void mm_canvas(std::string name = "v2_all.root") {
  TF1 *fitFcn = new TF1("fitFcn", fitFunction, 0.8, 1.2, 6);

  fitFcn->SetNpx(200);
  fitFcn->SetParNames("b1", "b2", "b3");
  fitFcn->SetParameters(1, 1, 1, 6, .03, 1);
  fitFcn->Update();

  TFile *f = new TFile(name.c_str());

  TCanvas *c1 = new TCanvas("c1", "c1", 1600, 900);

  TH2D *h1 = (TH2D *)f->Get("Missing_Mass/Missing_Mass");
  h1->Fit(fitFcn, "QM+");

  double l = 0.8;
  double h = 1.1;

  TLine *lower = new TLine(l, 0, l, 180000);
  TLine *upper = new TLine(h, 0, h, 180000);
  lower->SetLineColor(46);
  upper->SetLineColor(46);

  lower->SetLineWidth(2);
  upper->SetLineWidth(2);

  h1->Draw();
  lower->Draw("same");
  upper->Draw("same");
}

#if not defined(__CLING__)
int main(int argc, char const *argv[]) {
  if (argc < 2) {
    std::cerr << "Not enough arguments" << std::endl;
    std::cerr << "To Use:\t" << argv[0] << " v2_all.root" << std::endl;
    exit(1);
  }

  mm_canvas(argv[1]);

  return 0;
}
#endif