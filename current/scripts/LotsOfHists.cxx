#include <iostream>
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "THn.h"
#include "THnSparse.h"
#include "TStyle.h"

void LotsOfHists(const std::string &data_root, const std::string &mc_root) {
  TFile *root_data = new TFile(data_root.c_str());
  TFile *root_mc = new TFile(mc_root.c_str());
  TFile *out = new TFile("output.root", "RECREATE");
  out->cd();

  THnSparse *ndHist = (THnSparse *)root_data->Get("ndhist");
  // THnSparse *ndHist_rec = (THnSparse *)root_mc->Get("ndhist");
  // THnSparse *ndHist_thrown = (THnSparse *)root_mc->Get("ndhist_mc");
  // ndHist_rec->Divide(ndHist_thrown);
  // ndHist->Multiply(ndHist_rec);

  const int DIMENSIONS = ndHist->GetNdimensions();
  int nbins[DIMENSIONS];
  double xmin[DIMENSIONS];
  double xmax[DIMENSIONS];
  int bin[DIMENSIONS];
  for (int i = 0; i < DIMENSIONS; i++) nbins[i] = ndHist->GetAxis(i)->GetNbins();
  for (int i = 0; i < DIMENSIONS; i++) xmin[i] = ndHist->GetAxis(i)->GetBinLowEdge(1);
  for (int i = 0; i < DIMENSIONS; i++) xmax[i] = ndHist->GetAxis(i)->GetBinUpEdge(ndHist->GetAxis(i)->GetNbins());

  ndHist->Sumw2();
  TH1D *All_hists[nbins[0]][nbins[1]][nbins[0]];
  for (int w = 0; w < nbins[0]; w++) {
    for (int q2 = 0; q2 < nbins[1]; q2++) {
      for (int theta = 0; theta < nbins[2]; theta++) {
        bin[0] = w;
        bin[1] = q2;
        bin[2] = theta;
        float Q2_val = q2 * ((xmax[1] - xmin[1]) / nbins[1] * 1.0) + xmin[1];
        float W_val = w * ((xmax[0] - xmin[0]) / nbins[0] * 1.0) + xmin[0];
        float Theta_val = theta * ((xmax[2] - xmin[2]) / nbins[2] * 1.0) + xmin[2];

        All_hists[w][q2][theta] = new TH1D(Form("w%f_q2%f_theta%f", W_val, Q2_val, Theta_val),
                                           Form("w%f_q2%f_theta%f", W_val, Q2_val, Theta_val), nbins[3], -360, 360);
        for (int phi = 0; phi < nbins[3]; phi++) {
          bin[3] = phi;
          if (ndHist->GetBinContent(bin) == 0) continue;
          All_hists[w][q2][theta]->SetBinContent(phi, static_cast<double>(ndHist->GetBinContent(bin)));
          All_hists[w][q2][theta]->SetBinError(phi, static_cast<double>(ndHist->GetBinError(ndHist->GetBin(bin))));
        }
        if (All_hists[w][q2][theta]->GetEntries() > 4) {
          All_hists[w][q2][theta]->Fit("pol4");
          All_hists[w][q2][theta]->Write();
        }
      }
    }
  }
  out->Write();
  out->Close();
}

#if not defined(__CLING__)
int main(int argc, char const *argv[]) {
  if (argc < 2) {
    std::cerr << "Not enough arguments" << std::endl;
    std::cerr << "To Use:\t" << argv[0] << " data.root mc.root" << std::endl;
    exit(1);
  }

  LotsOfHists(argv[1], argv[2]);

  return 0;
}
#endif