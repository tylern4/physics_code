/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#ifndef HISTOGRAM_H
#define HISTOGRAM_H
#include <cmath>
#include <fstream>
#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THn.h"
#include "color.hpp"
#include "fits.hpp"
#include "missing_mass.hpp"

class mcHistogram {
 private:
 public:
  mcHistogram(std::string output_file);
  ~mcHistogram();
  TFile *RootOutputFile;
  TCanvas *def;
  int bins = 500;
  double p_min = 0.0;
  double p_max = 5.0;
  char hname[50];
  char htitle[500];

  static const int sector = 6;

  // Missing Mass
  int bins_MM = 300;
  double MM_min = 0.0;
  double MM_max = 3.0;
  // Missing Mass

  // W and Q^2
  double w_min = 0;
  double w_max = 3.25;
  double q2_min = 0;
  double q2_max = 5;

  static constexpr int W_bins = 40;
  static constexpr int Q2_bins = 5;
  static constexpr int theta_bins = 100;
  static constexpr int phi_bins = 100;
  double w_binned_min = 1.0;
  double w_binned_max = 2.0;
  double q2_binned_min = 1.0;
  double q2_binned_max = 4.0;

  double W_width = (w_binned_max - w_binned_min) / (double)W_bins;
  double Q2_width = (q2_binned_max - q2_binned_min) / (double)Q2_bins;

  TH2D *WvsQ2_hist = new TH2D("WvsQ2_hist", "W vs Q^{2}", bins, w_min, w_max, bins, q2_min, q2_max);
  TH1D *W_hist = new TH1D("W", "W", bins, w_min, w_max);

  TH2D *WvsQ2_MC = new TH2D("WvsQ2_MC", "W vs Q^{2} #pi^{+} N", bins, w_min, w_max, bins, q2_min, q2_max);
  TH1D *W_MC = new TH1D("W_MC", "W #pi^{+} N", bins, w_min, w_max);

  TH2D *WvsQ2_binned = new TH2D("WvsQ2_hist_binned", "W vs Q^{2} binned", W_bins, w_binned_min, w_binned_max, Q2_bins,
                                q2_binned_min, q2_binned_max);
  TH2D *WvsQ2_binned_MC = new TH2D("WvsQ2_hist_binned_MC", "W vs Q^{2} binned", W_bins, w_binned_min, w_binned_max,
                                   Q2_bins, q2_binned_min, q2_binned_max);

  TH1D *W_binned[Q2_bins];
  TH1D *W_binned_MC[Q2_bins];
  TH1D *Missing_Mass = new TH1D("Missing_Mass", "Missing Mass", bins_MM, MM_min, MM_max);
  TH1D *Missing_Mass_square = new TH1D("Missing_Mass_square", "Missing Mass^2", bins_MM, MM_min, MM_max);

  // W and Q^2
  void makeHists_W();
  void Fill_WQ2(double W, double Q2);
  void Fill_WQ2_MC(double W, double Q2);
  void WvsQ2_Write();
  void WvsQ2_binned_Write();

  // Missing Mass
  void Fill_Missing_Mass(MissingMass *miss_mass);
  void Write_Missing_Mass();
};

#endif
