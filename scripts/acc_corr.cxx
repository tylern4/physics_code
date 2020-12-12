#include <iostream>
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "THn.h"
#include "THnSparse.h"
#include "TStyle.h"

//////////////// W , Q2, Theta_star_pip, Phi_star_pip
void acc_corr(const std::string &data_root, const std::string &mc_root) {
  // Open root files
  TFile *root_data = new TFile(data_root.c_str());
  TFile *root_mc = new TFile(mc_root.c_str());

  // Get THnD from each file
  THnSparse *ndHist = (THnSparse *)root_data->Get("ndhist");
  THnSparse *ndHist_rec = (THnSparse *)root_mc->Get("ndhist");
  THnSparse *ndHist_thrown = (THnSparse *)root_mc->Get("ndhist_mc");
  // calculate acceptance
  ndHist_rec->Divide(ndHist_thrown);
  ndHist->Multiply(ndHist_rec);
  // Now ndHist shouild be acceptance corrected yeilds

  // Get number of dimentions from THnD
  const int DIMENSIONS = ndHist->GetNdimensions();
  // Make arrays of bins the length of the dimentions of the THnD
  // Since all the dimentions have their own number of bins assocatied with it
  // i.e. W has 25 bins, Q2 has 10 bins etc...
  int nbins[DIMENSIONS];

  // Make places to store the min and max values of the edge of the bins
  float xmin[DIMENSIONS];
  float xmax[DIMENSIONS];

  // Loop over the number of dimentions and
  for (int i = 0; i < DIMENSIONS; i++) {
    nbins[i] = ndHist->GetAxis(i)->GetNbins();
    xmin[i] = ndHist->GetAxis(i)->GetBinLowEdge(1);
    xmax[i] = ndHist->GetAxis(i)->GetBinUpEdge(ndHist->GetAxis(i)->GetNbins());
  }

  std::cout << "q2,w,theta,phi,y,yerr" << std::endl;
  int bin[DIMENSIONS];
  // Loop of number of bins in w,q2,theta,phi
  for (int w = 0; w < nbins[0]; w++) {
    for (int q2 = 0; q2 < nbins[1]; q2++) {
      for (int theta = 0; theta < nbins[2]; theta++) {
        for (int phi = 0; phi < nbins[3]; phi++) {
          // Set the values of the bin array to the values we want to extract
          bin[0] = w;
          bin[1] = q2;
          bin[2] = theta;
          bin[3] = phi;
          // If there is nothing in the bin just skip it
          if (ndHist->GetBinContent(bin) == 0) continue;
          // calculate and print out the bin centers for w,q2,theta,phi,yeilds,error
          std::cout << q2 * ((xmax[1] - xmin[1]) / nbins[1] * 1.0) + xmin[1] << ",";
          std::cout << w * ((xmax[0] - xmin[0]) / nbins[0] * 1.0) + xmin[0] << ",";
          std::cout << theta * ((xmax[2] - xmin[2]) / nbins[2] * 1.0) + xmin[2] << ",";
          std::cout << phi * ((xmax[3] - xmin[3]) / nbins[3] * 1.0) + xmin[3] << ",";
          std::cout << static_cast<float>(ndHist->GetBinContent(bin)) << ",";
          std::cout << static_cast<float>(ndHist->GetBinError(ndHist->GetBin(bin))) << std::endl;
        }
      }
    }
  }
}

// Hack to compile if not running at `root -l script.cxx`
#if not defined(__CLING__)
int main(int argc, char const *argv[]) {
  if (argc < 2) {
    std::cerr << "Not enough arguments" << std::endl;
    std::cerr << "To Use:\t" << argv[0] << " data.root mc.root" << std::endl;
    exit(1);
  }

  acc_corr(argv[1], argv[2]);

  return 0;
}
#endif