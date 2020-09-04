#include <iostream>
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"

TCanvas* EC_show(std::string infile) {
  TCanvas* c1 = new TCanvas("EC");

  TFile in_file(infile.c_str());
  TH2D* h;
  in_file.GetObject("EC_hists/EC_sampling_fraction_cut", h);

  TF1* fit1;
  in_file.GetObject("EC_hists/EC_P_fit", fit1);

  TF1* fit2;
  in_file.GetObject("EC_hists/EC_M_fit", fit2);

  TGraph* P;
  in_file.GetObject("EC_hists/Positive_EC_graph", P);

  TGraph* M;
  in_file.GetObject("EC_hists/Negative_EC_graph", M);

  h->Draw("colz");

  // fit1->Draw("same");
  // fit2->Draw("same");
  // P->Draw("*same");
  // M->Draw("*same");

  return c1;
}

#if !defined(__CLING__)
int main(int argc, char const* argv[]) {
  if (argc < 2) {
    std::cerr << "Not enough arguments" << std::endl;
    std::cerr << "To Use:\tEC_show data.root" << std::endl;
    exit(1);
  }

  auto can = EC_show(argv[1]);
  can->SaveAs("EC_show.pdf");
  return 0;
}
#endif