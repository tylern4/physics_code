#include <iostream>
#include "Math/Minimizer.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "THn.h"
#include "THnSparse.h"
#include "TStyle.h"

double TwoPhi(Double_t *x, Double_t *par) {
  float _phi = x[0] * M_PI / 180.0;
  float _epsilon = par[0];
  float _sigma_t = par[1];
  float _sigma_l = par[2];
  float _sigma_tt = par[3];
  float _sigma_lt = par[4];
  float ret = 0;

  ret = _sigma_t + _epsilon * _sigma_l;
  ret += _epsilon * _sigma_tt * cosf(2 * _phi);
  ret += sqrtf(2 * _epsilon * (1 + _epsilon)) * _sigma_lt * cosf(_phi);

  return static_cast<double>(ret);
}

void LotsOfHists(const std::string &data_root, const std::string &mc_root) {
  std::cout << "q2,w,theta,phi,y,yerr" << std::endl;
  TFile *root_data = new TFile(data_root.c_str());
  TFile *root_mc = new TFile(mc_root.c_str());
  TFile *out = new TFile("LotsOfHists.root", "RECREATE");
  out->cd();

  THnSparse *ndHist = (THnSparse *)root_data->Get("ndhist");
  THnSparse *ndHist_rec = (THnSparse *)root_mc->Get("ndhist");
  THnSparse *ndHist_thrown = (THnSparse *)root_mc->Get("ndhist_mc");
  ndHist_rec->Divide(ndHist_thrown);
  ndHist->Multiply(ndHist_rec);

  const int DIMENSIONS = ndHist->GetNdimensions();
  int nbins[DIMENSIONS];
  double xmin[DIMENSIONS];
  double xmax[DIMENSIONS];
  int bin[DIMENSIONS];
  for (int i = 0; i < DIMENSIONS; i++) nbins[i] = ndHist->GetAxis(i)->GetNbins();
  for (int i = 0; i < DIMENSIONS; i++) xmin[i] = ndHist->GetAxis(i)->GetBinLowEdge(1);
  for (int i = 0; i < DIMENSIONS; i++) xmax[i] = ndHist->GetAxis(i)->GetBinUpEdge(ndHist->GetAxis(i)->GetNbins());

  // ndHist->Sumw2();
  TH1D *All_hists[nbins[0]][nbins[1]][nbins[0]];
  TF1 *func = new TF1("func", TwoPhi, -360, 360, 5);
  func->SetParNames("#epsilon", "#sigma_{l}", "#sigma_{t}", "#sigma_{lt}", "#sigma_{tt}");
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  for (int w = 0; w < nbins[0]; w++) {
    for (int q2 = 0; q2 < nbins[1]; q2++) {
      for (int theta = 0; theta < nbins[2]; theta++) {
        bin[0] = w;
        bin[1] = q2;
        bin[2] = theta;
        float Q2_val = q2 * ((xmax[1] - xmin[1]) / nbins[1] * 1.0) + xmin[1];
        float W_val = w * ((xmax[0] - xmin[0]) / nbins[0] * 1.0) + xmin[0];
        float Theta_val = theta * ((xmax[2] - xmin[2]) / nbins[2] * 1.0) + xmin[2];

        All_hists[w][q2][theta] =
            new TH1D(Form("w%0.2f_qSq%0.2f_theta%0.2f", W_val, Q2_val, Theta_val),
                     Form("w%0.2f_qSq%0.2f_theta%0.2f", W_val, Q2_val, Theta_val), nbins[3], -360, 360);
        for (int phi = 0; phi < nbins[3]; phi++) {
          bin[3] = phi;
          float Phi_val = phi * ((xmax[3] - xmin[3]) / nbins[3] * 1.0) + xmin[3];
          if (ndHist->GetBinContent(bin) == 0) continue;
          std::cout << W_val << "," << Q2_val << "," << Theta_val << "," << Phi_val << ","
                    << static_cast<double>(ndHist->GetBinContent(bin)) << ","
                    << static_cast<double>(ndHist->GetBinError(ndHist->GetBin(bin))) << std::endl;
          All_hists[w][q2][theta]->SetBinContent(phi, static_cast<double>(ndHist->GetBinContent(bin)));
          All_hists[w][q2][theta]->SetBinError(phi, static_cast<double>(ndHist->GetBinError(ndHist->GetBin(bin))));
        }
        if (All_hists[w][q2][theta]->GetEntries() > 4) {
          // std::cout << w << "\t" << q2 << "\t" << theta << std::endl;
          // double par[5] = {1.0, 1.0, 1.0, 1.0, 1.0};
          // func->SetParameters(par);
          // for (int i = 0; i < 10; i++) All_hists[w][q2][theta]->Fit("func", "QMN");
          // All_hists[w][q2][theta]->Fit("func", "QM+");
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