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

// a + b cos(phi) + c cos(2phi)
double TwoPhi(Double_t *x, Double_t *par) {
  double _phi = x[0];  // * M_PI / 180.0;
  double _a = par[0];
  double _b = par[1];
  double _c = par[2];

  double ret = 0;
  ret = _a + _b * cos(_phi) + _c * cos(2.0 * _phi);

  return static_cast<double>(ret);
}

void LotsOfHists(const std::string &data_root, const std::string &mc_root) {
  std::cout << "w,q2,cos_theta,phi,y,yerr" << std::endl;
  TFile *root_data = new TFile(data_root.c_str());
  TFile *root_mc = new TFile(mc_root.c_str());
  TFile *out = new TFile("LotsOfHists.root", "RECREATE");
  out->cd();

  THnD *ndHist = (THnD *)root_data->Get("ndhist");
  THnD *ndHist_rec = (THnD *)root_mc->Get("ndhist");
  THnD *ndHist_acc = (THnD *)root_mc->Get("ndhist");
  THnD *ndHist_thrown = (THnD *)root_mc->Get("ndhist_mc");

  ndHist->Sumw2();
  ndHist_rec->Sumw2();
  ndHist_thrown->Sumw2();

  const int DIMENSIONS = ndHist->GetNdimensions();
  for (size_t i = 0; i < DIMENSIONS; i++) {
    auto single_acc = ndHist_rec->Projection(i);
    single_acc->Divide(ndHist_thrown->Projection(i));
    single_acc->Multiply(ndHist->Projection(i));
    single_acc->Write(Form("Acceptance_%zu", i));
  }

  // ndHist_thrown->Divide(ndHist_rec);
  ndHist_rec->Divide(ndHist_thrown);
  ndHist_rec->Write("Acceptance");
  // ndHist_acc->Divide(ndHist_thrown);
  // ndHist_acc->Write("Acceptance");
  // plot acceptane histograms as well
  ndHist->Multiply(ndHist_rec);
  // ndHist->Divide(ndHist_rec);

  int nbins[DIMENSIONS];
  double xmin[DIMENSIONS];
  double xmax[DIMENSIONS];
  int bin[DIMENSIONS];
  for (int i = 0; i < DIMENSIONS; i++) nbins[i] = ndHist->GetAxis(i)->GetNbins();
  for (int i = 0; i < DIMENSIONS; i++) xmin[i] = ndHist->GetAxis(i)->GetBinLowEdge(1);
  for (int i = 0; i < DIMENSIONS; i++) xmax[i] = ndHist->GetAxis(i)->GetBinUpEdge(ndHist->GetAxis(i)->GetNbins());

  ndHist->Sumw2();
  TH1D *All_hists[nbins[0]][nbins[1]][nbins[2]];
  TF1 *func = new TF1("func", TwoPhi, -360, 360, 3);
  // func->SetParNames("#epsilon", "#sigma_{l}", "#sigma_{t}", "#sigma_{lt}", "#sigma_{tt}");
  func->SetParNames("#alpha", "#beta", "#gamma");
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
            new TH1D(Form("w%0.5f_qSq%0.5f_theta%0.5f", W_val, Q2_val, Theta_val),
                     Form("w%0.5f_qSq%0.5f_theta%0.5f", W_val, Q2_val, Theta_val), nbins[3], xmin[3], xmax[3]);
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

        if (All_hists[w][q2][theta]->GetEntries() >= 5) {
          // std::cout << w << "\t" << q2 << "\t" << theta << std::endl;
          double par[3] = {1.0, 1.0, 1.0};
          func->SetParameters(par);
          for (int i = 0; i < 10; i++) All_hists[w][q2][theta]->Fit("func", "QMN");
          All_hists[w][q2][theta]->Fit("func", "QM+");
          // All_hists[w][q2][theta]->Write();
        } else {
          delete All_hists[w][q2][theta];
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