/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

// Reference: https://pprc.qmul.ac.uk/~bevan/yeti/fitting.pdf

#include "fits.hpp"

void Fits::FitGaus(TH1D *hist, double min_value, double max_value) {
  if (hist->GetEntries() > 1000) {
    // ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    TF1 *fitFunc = new TF1("fitFunc", gaus.c_str(), min_value, max_value);
    fitFunc->SetLineColor(2);
    par_max = std::isnan(hist->GetMaximum()) ? 0 : hist->GetMaximum();
    par_mean = std::isnan(hist->GetMean()) ? 0 : hist->GetMean();
    fitFunc->SetParameter(0, par_max);
    fitFunc->SetParameter(1, par_mean);
    fitFunc->SetParameter(2, 1);
    fitFunc->SetParNames("height", "mean", "FWHM");

    hist->Fit("fitFunc", "q", "", min_value, max_value);

    par_mean = std::isnan(fitFunc->GetParameter("mean"))
                   ? 0
                   : fitFunc->GetParameter("mean");
    par_FWHM = std::isnan(fitFunc->GetParameter("FWHM"))
                   ? 0
                   : fitFunc->GetParameter("FWHM");

    fitFunc->SetParameter(0, par_max);
    fitFunc->SetParameter(1, par_mean);
    fitFunc->SetParameter(2, par_FWHM);
    for (int i = 0; i < 10; i++)
      hist->Fit("fitFunc", "q", "", min_value, max_value);

    mean = fitFunc->GetParameter("mean");
    FWHM = fitFunc->GetParameter("FWHM");
    sigma = fitFunc->GetParameter("FWHM") /
            (2 * sqrt(2 * log(2))); // 2.35482004503;
  }
}

void Fits::FitLandau(TH1D *hist, double min_value, double max_value) {
  if (hist->GetEntries() > 1000) {
    // ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    TF1 *fitFunc = new TF1("fitFunc", "landau", min_value, max_value);

    hist->Fit("fitFunc", "q", "", min_value, max_value);

    for (int i = 0; i < 10; i++)
      hist->Fit("fitFunc", "q", "", min_value, max_value);
  }
}

void Fits::FitPoly_1D(TH1D *hist, double min_value, double max_value) {
  // ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  TF1 *fitFunc = new TF1("fitFunc", ploy1d.c_str(), min_value, max_value);
  fitFunc->SetLineColor(9);
  fitFunc->SetParNames("a", "b");

  hist->Fit("fitFunc", "q", "", min_value, max_value);

  fitFunc->SetParameter(0, fitFunc->GetParameter("a"));
  fitFunc->SetParameter(1, fitFunc->GetParameter("b"));

  hist->Fit("fitFunc", "q", "", min_value, max_value);

  a = fitFunc->GetParameter("a");
  b = fitFunc->GetParameter("b");
}

void Fits::FitPoly_2D(TH1D *hist, double min_value, double max_value) {
  // ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  TF1 *fitFunc = new TF1("fitFunc", ploy2d.c_str(), min_value, max_value);
  fitFunc->SetLineColor(30);
  fitFunc->SetParNames("a", "b", "c");

  hist->Fit("fitFunc", "q", "", min_value, max_value);

  fitFunc->SetParameter(0, fitFunc->GetParameter("a"));
  fitFunc->SetParameter(1, fitFunc->GetParameter("b"));
  fitFunc->SetParameter(2, fitFunc->GetParameter("c"));

  hist->Fit("fitFunc", "q", "", min_value, max_value);

  a = fitFunc->GetParameter("a");
  b = fitFunc->GetParameter("b");
  c = fitFunc->GetParameter("c");
}

void Fits::FitPoly_3D(TH1D *hist, double min_value, double max_value) {
  // ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  TF1 *fitFunc = new TF1("fitFunc", ploy3d.c_str(), min_value, max_value);
  fitFunc->SetLineColor(46);
  fitFunc->SetParNames("a", "b", "c", "d");

  hist->Fit("fitFunc", "q", "", min_value, max_value);

  fitFunc->SetParameter(0, fitFunc->GetParameter("a"));
  fitFunc->SetParameter(1, fitFunc->GetParameter("b"));
  fitFunc->SetParameter(2, fitFunc->GetParameter("c"));
  fitFunc->SetParameter(3, fitFunc->GetParameter("d"));

  hist->Fit("fitFunc", "q", "", min_value, max_value);

  a = fitFunc->GetParameter("a");
  b = fitFunc->GetParameter("b");
  c = fitFunc->GetParameter("c");
  d = fitFunc->GetParameter("d");
}

void Fits::FitPoly_4D(TH1D *hist, double min_value, double max_value) {
  // ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  TF1 *fitFunc = new TF1("fitFunc", ploy4d.c_str(), min_value, max_value);
  fitFunc->SetLineColor(42);
  fitFunc->SetParNames("a", "b", "c", "d", "e");

  hist->Fit("fitFunc", "q", "", min_value, max_value);

  fitFunc->SetParameter(0, fitFunc->GetParameter("a"));
  fitFunc->SetParameter(1, fitFunc->GetParameter("b"));
  fitFunc->SetParameter(2, fitFunc->GetParameter("c"));
  fitFunc->SetParameter(3, fitFunc->GetParameter("d"));
  fitFunc->SetParameter(4, fitFunc->GetParameter("e"));

  hist->Fit("fitFunc", "q", "", min_value, max_value);

  a = fitFunc->GetParameter("a");
  b = fitFunc->GetParameter("b");
  c = fitFunc->GetParameter("c");
  d = fitFunc->GetParameter("d");
  e = fitFunc->GetParameter("e");
}

void Fits::FitPoly_fid(TH2D *hist, double min_value, double max_value) {
  // ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  TF1 *fitFunc = new TF1("fitFunc", "pol8", min_value, max_value);
  fitFunc->SetLineColor(42);
  fitFunc->SetParNames("a", "b", "c", "d", "e", "f", "g", "h");

  hist->Fit("fitFunc", "q", "", min_value, max_value);

  fitFunc->SetParameter(0, fitFunc->GetParameter("a"));
  fitFunc->SetParameter(1, fitFunc->GetParameter("b"));
  fitFunc->SetParameter(2, fitFunc->GetParameter("c"));
  fitFunc->SetParameter(3, fitFunc->GetParameter("d"));
  fitFunc->SetParameter(4, fitFunc->GetParameter("e"));
  fitFunc->SetParameter(5, fitFunc->GetParameter("f"));
  fitFunc->SetParameter(6, fitFunc->GetParameter("g"));
  fitFunc->SetParameter(7, fitFunc->GetParameter("h"));

  hist->Fit("fitFunc", "q", "", min_value, max_value);

  // a = fitFunc->GetParameter("a");
  // b = fitFunc->GetParameter("b");
  // c = fitFunc->GetParameter("c");
  // d = fitFunc->GetParameter("d");
  // e = fitFunc->GetParameter("e");
}

double Fits::fiducial_phi_lo(double theta_e, double theta_e_min, double k,
                             double m, int c) {
  return -c * pow(sin((theta_e - theta_e_min) * 0.01745),
                  k + m / theta_e + 1500. / (theta_e * theta_e));
}

double Fits::fiducial_phi_hi(double theta_e, double theta_e_min, double k,
                             double m, int c) {
  return c * pow(sin((theta_e - theta_e_min) * 0.01745),
                 k + m / theta_e + 1500. / (theta_e * theta_e));
}

void Fits::FitFiducial_lo(TH2D *hist2d, double min_value, double max_value) {
  // ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  a = b = c = d = 0.5;
  TF1 *fitFunc_lo =
      new TF1("fitFunc_lo", "-[0]*TMath::Power(TMath::Sin((x-[1])*0.01745),"
                            "[2] +[3] / x + 1500. / (x * x))",
              min_value, max_value);

  fitFunc_lo->SetLineColor(7);
  fitFunc_lo->SetParNames("a", "b", "c", "d");

  hist2d->Fit("fitFunc_lo", "q", "", min_value, max_value);

  fitFunc_lo->SetParameter(0, fitFunc_lo->GetParameter("a"));
  fitFunc_lo->SetParameter(1, fitFunc_lo->GetParameter("b"));
  fitFunc_lo->SetParameter(2, fitFunc_lo->GetParameter("c"));
  fitFunc_lo->SetParameter(3, fitFunc_lo->GetParameter("d"));

  hist2d->Fit("fitFunc_lo", "q", "", min_value, max_value);

  a = fitFunc_lo->GetParameter("a");
  b = fitFunc_lo->GetParameter("b");
  c = fitFunc_lo->GetParameter("c");
  d = fitFunc_lo->GetParameter("d");
}

void Fits::FitFiducial_hi(TH2D *hist2d, double min_value, double max_value) {
  // ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  a = b = c = d = 0.5;
  TF1 *fitFunc_hi =
      new TF1("fitFunc_hi", "[0]*TMath::Power(TMath::Sin((x-[1])*"
                            "0.01745),[2] +[3] / x + 1500. / (x * x))",
              max_value, max_value);

  fitFunc_hi->SetLineColor(8);
  fitFunc_hi->SetParNames("a", "b", "c", "d");

  hist2d->Fit("fitFunc_hi", "q", "", max_value, max_value);

  fitFunc_hi->SetParameter(0, fitFunc_hi->GetParameter("a"));
  fitFunc_hi->SetParameter(1, fitFunc_hi->GetParameter("b"));
  fitFunc_hi->SetParameter(2, fitFunc_hi->GetParameter("c"));
  fitFunc_hi->SetParameter(3, fitFunc_hi->GetParameter("d"));

  hist2d->Fit("fitFunc_hi", "q", "", max_value, max_value);

  a = fitFunc_hi->GetParameter("a");
  b = fitFunc_hi->GetParameter("b");
  c = fitFunc_hi->GetParameter("c");
  d = fitFunc_hi->GetParameter("d");
}

void Fits::FitFiducial(TH2D *hist2d, double min_value, double max_value) {
  // ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  a = b = c = d = 0.5;
  TF1 *fitFunc = new TF1("fitFunc", "[0]*TMath::Power(TMath::Sin((x-[1])*"
                                    "0.01745),[2] +[3] / x + 1500. / (x * x))",
                         max_value, max_value);

  fitFunc->SetLineColor(9);
  fitFunc->SetParNames("a", "b", "c", "d");

  hist2d->Fit("fitFunc", "q", "", max_value, max_value);

  fitFunc->SetParameter(0, fitFunc->GetParameter("a"));
  fitFunc->SetParameter(1, fitFunc->GetParameter("b"));
  fitFunc->SetParameter(2, fitFunc->GetParameter("c"));
  fitFunc->SetParameter(3, fitFunc->GetParameter("d"));

  hist2d->Fit("fitFunc", "q", "", max_value, max_value);

  a = fitFunc->GetParameter("a");
  b = fitFunc->GetParameter("b");
  c = fitFunc->GetParameter("c");
  d = fitFunc->GetParameter("d");
}

void Fits::FitGenNormal(TH1D *hist, double min_value, double max_value) {
  double min, max, val, min_m, max_m;
  if (hist->GetEntries() > 1000) {
    // ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit");

    TF1 *fitFunc = new TF1("genNormal", genNormal, min_value, max_value, 4);

    fitFunc->SetParLimits(1, 5.0, 200.0);

    fitFunc->SetParameter(0, 15.0);
    fitFunc->SetParameter(1, 10.0);
    fitFunc->SetParameter(2, (min_value + max_value) / 2.0);
    fitFunc->SetParameter(3, 3000.0);
    fitFunc->SetParNames("alpha", "beta", "mu", "weight");

    hist->Fit("genNormal", "Q+", "", min_value, max_value);

    for (double m = min_value; m < max_value; m = m + 0.005) {
      val = fitFunc->Derivative(m);

      if (max < val) {
        max = val;
        max_m = m;
      }
      if (val < min) {
        min = val;
        min_m = m;
      }
    }
    min_edge_x = min_m;
    max_edge_x = max_m;
  }
}
