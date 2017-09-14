/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#include "fits.hpp"

void Fits::FitGaus(TH1D *hist, double min_value, double max_value) {
  TF1 *fitFunc = new TF1("fitFunc", gaus.c_str(), min_value, max_value);
  fitFunc->SetLineColor(2);
  par_max = std::isnan(hist->GetMaximum()) ? 0 : hist->GetMaximum();
  par_mean = std::isnan(hist->GetMean()) ? 0 : hist->GetMean();
  fitFunc->SetParameter(0, par_max);
  fitFunc->SetParameter(1, par_mean);
  fitFunc->SetParameter(2, 1);
  fitFunc->SetParNames("height", "mean", "FWHM");

  hist->Fit("fitFunc", "qM0+", "", min_value, max_value);

  par_mean = std::isnan(fitFunc->GetParameter("mean"))
                 ? 0
                 : fitFunc->GetParameter("mean");
  par_FWHM = std::isnan(fitFunc->GetParameter("FWHM"))
                 ? 0
                 : fitFunc->GetParameter("FWHM");

  fitFunc->SetParameter(0, par_max);
  fitFunc->SetParameter(1, par_mean);
  fitFunc->SetParameter(2, par_FWHM);
  hist->Fit("fitFunc", "qM+", "", min_value, max_value);

  mean = fitFunc->GetParameter("mean");
  FWHM = fitFunc->GetParameter("FWHM");
  sigma =
      fitFunc->GetParameter("FWHM") / (2 * sqrt(2 * log(2))); // 2.35482004503;
  gStyle->SetOptFit(1111);
}

void Fits::FitPoly_1D(TH1D *hist, double min_value, double max_value) {
  TF1 *fitFunc = new TF1("fitFunc", ploy1d.c_str(), min_value, max_value);
  fitFunc->SetLineColor(9);
  fitFunc->SetParNames("a", "b");

  hist->Fit("fitFunc", "qM0+", "", min_value, max_value);

  fitFunc->SetParameter(0, fitFunc->GetParameter("a"));
  fitFunc->SetParameter(1, fitFunc->GetParameter("b"));

  hist->Fit("fitFunc", "qM+", "", min_value, max_value);

  a = fitFunc->GetParameter("a");
  b = fitFunc->GetParameter("b");
  gStyle->SetOptFit(1111);
}

void Fits::FitPoly_2D(TH1D *hist, double min_value, double max_value) {
  TF1 *fitFunc = new TF1("fitFunc", ploy2d.c_str(), min_value, max_value);
  fitFunc->SetLineColor(30);
  fitFunc->SetParNames("a", "b", "c");

  hist->Fit("fitFunc", "qM0+", "", min_value, max_value);

  fitFunc->SetParameter(0, fitFunc->GetParameter("a"));
  fitFunc->SetParameter(1, fitFunc->GetParameter("b"));
  fitFunc->SetParameter(2, fitFunc->GetParameter("c"));

  hist->Fit("fitFunc", "qM+", "", min_value, max_value);

  a = fitFunc->GetParameter("a");
  b = fitFunc->GetParameter("b");
  c = fitFunc->GetParameter("c");
  gStyle->SetOptFit(1111);
}

void Fits::FitPoly_3D(TH1D *hist, double min_value, double max_value) {
  TF1 *fitFunc = new TF1("fitFunc", ploy3d.c_str(), min_value, max_value);
  fitFunc->SetLineColor(46);
  fitFunc->SetParNames("a", "b", "c", "d");

  hist->Fit("fitFunc", "qM0+", "", min_value, max_value);

  fitFunc->SetParameter(0, fitFunc->GetParameter("a"));
  fitFunc->SetParameter(1, fitFunc->GetParameter("b"));
  fitFunc->SetParameter(2, fitFunc->GetParameter("c"));
  fitFunc->SetParameter(3, fitFunc->GetParameter("d"));

  hist->Fit("fitFunc", "qM+", "", min_value, max_value);

  a = fitFunc->GetParameter("a");
  b = fitFunc->GetParameter("b");
  c = fitFunc->GetParameter("c");
  d = fitFunc->GetParameter("d");
  gStyle->SetOptFit(1111);
}

void Fits::FitPoly_4D(TH1D *hist, double min_value, double max_value) {
  TF1 *fitFunc = new TF1("fitFunc", ploy4d.c_str(), min_value, max_value);
  fitFunc->SetLineColor(42);
  fitFunc->SetParNames("a", "b", "c", "d", "e");

  hist->Fit("fitFunc", "qM0+", "", min_value, max_value);

  fitFunc->SetParameter(0, fitFunc->GetParameter("a"));
  fitFunc->SetParameter(1, fitFunc->GetParameter("b"));
  fitFunc->SetParameter(2, fitFunc->GetParameter("c"));
  fitFunc->SetParameter(3, fitFunc->GetParameter("d"));
  fitFunc->SetParameter(4, fitFunc->GetParameter("e"));

  hist->Fit("fitFunc", "qM+", "", min_value, max_value);

  a = fitFunc->GetParameter("a");
  b = fitFunc->GetParameter("b");
  c = fitFunc->GetParameter("c");
  d = fitFunc->GetParameter("d");
  e = fitFunc->GetParameter("e");
  gStyle->SetOptFit(1111);
}

double Fits::fiducial_phi_lo(double theta_e, double theta_e_min, double k,
                             double m) {
  return -37.14 * pow(sin((theta_e - theta_e_min) * 0.01745),
                      k + m / theta_e + 1500. / (theta_e * theta_e));
}

double Fits::fiducial_phi_hi(double theta_e, double theta_e_min, double k,
                             double m) {
  return 37.14 * pow(sin((theta_e - theta_e_min) * 0.01745),
                     k + m / theta_e + 1500. / (theta_e * theta_e));
}

void Fits::FitFiducial_lo(TH2D *hist2d, double min_value, double max_value) {
  a = b = c = d = 0.5;
  TF1 *fitFunc_lo =
      new TF1("fitFunc_lo", "-37.14*TMath::Power(TMath::Sin(([0]-[1])*0.01745),"
                            "[2] +[3] / [0] + 1500. / ([0] * [0]))",
              min_value, 1 + ((min_value + max_value) / 2.0));

  fitFunc_lo->SetLineColor(41);
  fitFunc_lo->SetParNames("a", "b", "c", "d");

  hist2d->Fit("fitFunc_lo", "qM0+", "", min_value,
              (min_value + max_value) / 2.0);

  fitFunc_lo->SetParameter(0, fitFunc_lo->GetParameter("a"));
  fitFunc_lo->SetParameter(1, fitFunc_lo->GetParameter("b"));
  fitFunc_lo->SetParameter(2, fitFunc_lo->GetParameter("c"));
  fitFunc_lo->SetParameter(3, fitFunc_lo->GetParameter("d"));
  std::cout << a << "\t" << b << "\t" << c << "\t" << d << "\t" << std::endl;
  hist2d->Fit("fitFunc_lo", "qM+", "", min_value,
              (min_value + max_value) / 2.0);

  a = fitFunc_lo->GetParameter("a");
  b = fitFunc_lo->GetParameter("b");
  c = fitFunc_lo->GetParameter("c");
  d = fitFunc_lo->GetParameter("d");
  std::cout << a << "\t" << b << "\t" << c << "\t" << d << "\t" << std::endl;
  gStyle->SetOptFit(1111);
}

void Fits::FitFiducial_hi(TH2D *hist2d, double min_value, double max_value) {
  a = b = c = d = 0.5;
  TF1 *fitFunc_hi =
      new TF1("fitFunc_hi", "37.14*TMath::Power(TMath::Sin(([0]-[1])*"
                            "0.01745),[2] +[3] / [0] + 1500. / ([0] * [0]))",
              ((min_value + max_value) / 2.0) - 1, max_value);

  fitFunc_hi->SetLineColor(42);
  fitFunc_hi->SetParNames("a", "b", "c", "d");

  hist2d->Fit("fitFunc_hi", "qM0+", "", (min_value + max_value) / 2.0,
              max_value);

  fitFunc_hi->SetParameter(0, fitFunc_hi->GetParameter("a"));
  fitFunc_hi->SetParameter(1, fitFunc_hi->GetParameter("b"));
  fitFunc_hi->SetParameter(2, fitFunc_hi->GetParameter("c"));
  fitFunc_hi->SetParameter(3, fitFunc_hi->GetParameter("d"));

  hist2d->Fit("fitFunc_hi", "qM+", "", (min_value + max_value) / 2.0,
              max_value);

  a = fitFunc_hi->GetParameter("a");
  b = fitFunc_hi->GetParameter("b");
  c = fitFunc_hi->GetParameter("c");
  d = fitFunc_hi->GetParameter("d");

  gStyle->SetOptFit(1111);
}
