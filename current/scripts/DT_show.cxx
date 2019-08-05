#include <iostream>
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"

TCanvas* DT_show(std::string infile = "v2_all.root") {
  TFile in_file(infile.c_str());
  TCanvas* c1 = new TCanvas();
  TH2D* Proton_h;
  in_file.GetObject("Delta_T/delta_t_mass_P", Proton_h);

  TF1* Proton_fit1;
  in_file.GetObject("Delta_T/Proton_Pos_fit", Proton_fit1);

  TF1* Proton_fit2;
  in_file.GetObject("Delta_T/Proton_Neg_fit", Proton_fit2);

  TGraph* Proton_graph1;
  in_file.GetObject("Delta_T/Proton_Pos_graph", Proton_graph1);

  TGraph* Proton_graph2;
  in_file.GetObject("Delta_T/Proton_Neg_graph", Proton_graph2);

  c1->cd();
  Proton_h->Draw();
  Proton_fit1->Draw("same");
  Proton_fit2->Draw("same");
  Proton_graph1->Draw("*same");
  Proton_graph2->Draw("*same");

  TCanvas* c2 = new TCanvas();
  TH2D* Pip_h;
  in_file.GetObject("Delta_T/delta_t_mass_PIP", Pip_h);

  TF1* Pip_fit1;
  in_file.GetObject("Delta_T/Pip_Pos_fit", Pip_fit1);

  TF1* Pip_fit2;
  in_file.GetObject("Delta_T/Pip_Neg_fit", Pip_fit2);

  TGraph* Pip_graph1;
  in_file.GetObject("Delta_T/Pip_Pos_graph", Pip_graph1);

  TGraph* Pip_graph2;
  in_file.GetObject("Delta_T/Pip_Neg_graph", Pip_graph2);

  c2->cd();
  Pip_h->Draw();
  Pip_fit1->Draw("same");
  Pip_fit2->Draw("same");
  Pip_graph1->Draw("Psame");
  Pip_graph2->Draw("Psame");

  return c1;
}

#if !defined(__CLING__)
int main(int argc, char const* argv[]) {
  if (argc < 2) {
    std::cerr << "Not enough arguments" << std::endl;
    std::cerr << "To Use:\tDT_show data.root" << std::endl;
    exit(1);
  }

  auto can = DT_show(argv[1]);
  can->SaveAs("DT_show.pdf");
  return 0;
}
#endif