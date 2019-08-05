#include <iostream>
#include "TCanvas.h"
#include "TFile.h"
#include "TH2D.h"
#include "THn.h"

static const int _W = 0;
static const int _Q2 = 1;
static const int _SECTOR = 2;
static const int _THETA = 3;
static const int _PHI = 4;

int theta_phi() {
  TFile *f = new TFile("today2.root");
  std::stringstream ss;
  TCanvas *c1 = new TCanvas("c1", "c1", 1600, 900);

  c1->Divide(5);
  double width = 0.6;
  for (int i = 0; i < 5; i++) {
    c1->cd(i + 1);
    std::stringstream ss;
    ss << "W_Q2_binned/W_" << std::fixed << std::setprecision(3) << (1.0 + i * 0.6) << "_" << (1.0 + (i + 1) * 0.6);
    std::string name = ss.str();
    f->Get(name.c_str())->Draw("");
  }
  return 0;
}
