/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef CUTS_HPP
#define CUTS_HPP
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"

class Cuts {
private:
  double par_max, par_mean, par_FWHM;
  std::string ploy1d = "[0]+[1]*x";
  std::string ploy2d = ploy1d + "+[2]*x*x";
  std::string ploy3d = ploy2d + "+[3]*x*x*x";
  std::string ploy4d = ploy3d + "+[4]*x*x*x*x";
  std::string gaus = "[0]*TMath::Gaus(x,[1],[2],1)";

public:
  double mean, sigma, FWHM;
  double a, b, c, d, e;
  void FitGaus(TH1D *hist, double min_value, double max_value) {
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
    sigma = fitFunc->GetParameter("FWHM") /
            (2 * sqrt(2 * log(2))); // 2.35482004503;
    gStyle->SetOptFit(1111);
  }

  void FitPoly_1D(TH1D *hist, double min_value, double max_value) {
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

  void FitPoly_2D(TH1D *hist, double min_value, double max_value) {
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

  void FitPoly_3D(TH1D *hist, double min_value, double max_value) {
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

  void FitPoly_4D(TH1D *hist, double min_value, double max_value) {
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
};

#endif