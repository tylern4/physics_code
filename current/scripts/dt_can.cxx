#include <iostream>
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"

TCanvas *dt_can(std::string infile) {
  TFile *f = new TFile(infile.c_str());

  TCanvas *c1 = new TCanvas("c1", "c1", 1600, 900);
  c1->SetLogz();
  TH2D *h1 = (TH2D *)f->Get("Delta_T/delta_t_mass_PIP");
  // gStyle->SetPalette(kCividis);
  // TColor::InvertPalette();
  h1->Draw();

  return c1;
}

#if !defined(__CLING__)
int main(int argc, char const *argv[]) {
  if (argc < 3) {
    std::cerr << "Not enough arguments" << std::endl;
    std::cerr << "To Use:\tdt_can data.root" << std::endl;
    exit(1);
  }

  auto can = dt_can(argv[1]);
  can->SaveAs("dt_can.pdf");
  return 0;
}
#endif