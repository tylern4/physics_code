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
#include <iostream>

void dt_can(
    std::string file_name = "/Users/tylern/Desktop/show/today_sim.root") {
  TFile *f = new TFile(file_name.c_str());

  TCanvas *c1 = new TCanvas("c1", "c1", 1600, 900);
  c1->SetLogz();
  TH2D *h1 = (TH2D *)f->Get("Delta_T/delta_t_mass_PIP");
  gStyle->SetPalette(kCividis);
  // TColor::InvertPalette();
  h1->GetXaxis()->SetRange(100, 200);
  h1->Draw();
}

#if not defined(__CLING__)
int main(int argc, char const *argv[]) {
  if (argc < 2) {
    std::cerr << "Not enough arguments" << std::endl;
    std::cerr << "To Use:\t" << argv[0] << " v2_all.root" << std::endl;
    exit(1);
  }

  dt_can(argv[1]);

  return 0;
}
#endif