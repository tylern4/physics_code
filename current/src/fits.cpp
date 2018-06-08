/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

// Reference: https://pprc.qmul.ac.uk/~bevan/yeti/fitting.pdf
#include "fits.hpp"

Fits::Fits() {}
Fits::~Fits() {}

void Fits::Set_min(double val) { min_value = val; };
void Fits::Set_max(double val) { max_value = val; };
void Fits::Set_lineColor(int val) { color = val; };

double Fits::Get_left_edge() { return left_edge_x; }
double Fits::Get_right_edge() { return right_edge_x; }
double Fits::Get_sigma() { return sigma; }
double Fits::Get_mean() { return mean; }
double Fits::Get_FWHM() { return FWHM; }

TF1 *Fits::FitGaus(TH1D *hist) {
  if (hist->GetEntries() > 100) {
    // if (hist->GetEntries() > 10000) ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    TF1 *fitFunc = new TF1("fitFunc", func::gausian, -100.0, 100.0, 3);
    // TF1 *fitFunc = new TF1("fitFunc", "gaus", min_value, max_value);
    fitFunc->SetLineColor(color);
    par_max = std::isnan(hist->GetMaximum()) ? 0 : hist->GetMaximum();
    par_mean = std::isnan(hist->GetMean()) ? 0 : hist->GetMean();
    par_RMS = std::isnan(hist->GetRMS()) ? 0 : hist->GetRMS();
    fitFunc->SetParameter(0, par_max);
    fitFunc->SetParameter(1, par_mean);
    fitFunc->SetParameter(2, par_RMS);
    fitFunc->SetParNames("constant", "mean", "#sigma");

    hist->Fit("fitFunc", "QM0+", "", min_value, max_value);
    for (size_t i = 0; i < 10; i++) {
      par_mean = std::isnan(fitFunc->GetParameter("mean")) ? 0 : fitFunc->GetParameter("mean");
      par_RMS = std::isnan(fitFunc->GetParameter("#sigma")) ? 0 : fitFunc->GetParameter("#sigma");
      fitFunc->SetParameter(0, par_max);
      fitFunc->SetParameter(1, par_mean);
      fitFunc->SetParameter(2, par_RMS);
      hist->Fit("fitFunc", "QM0+", "", min_value, max_value);
    }

    hist->Fit("fitFunc", "QM+", "", min_value, max_value);

    mean = fitFunc->GetParameter("mean");
    FWHM = fitFunc->GetParameter("#sigma");
    sigma = fitFunc->GetParameter("#sigma") / (2 * sqrt(2 * log(2)));  // 2.35482004503;

    left_edge_x = mean - 3 * sigma;
    right_edge_x = mean + 3 * sigma;

    return fitFunc;
  }
}

TF1 *Fits::FitLandauGaus(TH1D *hist) {
  if (hist->GetEntries() > 1000) {
    double par[6];
    // TF1 *fitFuncGaus = new TF1("fitFuncGaus", func::gausian, 30.0, 250.0, 3);
    TF1 *fitFuncGaus = new TF1("gaus", "gaus", 30.0, 250.0);
    TF1 *fitFuncLandau = new TF1("landau", "landau", 0.0, 30.0);
    TF1 *total = new TF1("total", "landau(0)+gaus(3)", 0.0, 250.0);

    fitFuncLandau->SetLineColor(9);
    fitFuncLandau->SetParameter(0, 1);
    fitFuncLandau->SetParameter(1, 1);
    fitFuncLandau->SetParameter(2, 1);
    hist->Fit("landau", "RQM+", "", 0.0, 30.0);

    fitFuncGaus->SetLineColor(8);
    par_max = std::isnan(hist->GetMaximum()) ? 0 : hist->GetMaximum();
    par_mean = std::isnan(hist->GetMean()) ? 0 : hist->GetMean();
    par_RMS = std::isnan(hist->GetRMS()) ? 0 : hist->GetRMS();
    fitFuncGaus->SetParameter(0, par_max);
    fitFuncGaus->SetParameter(1, par_mean);
    fitFuncGaus->SetParameter(2, par_RMS);
    fitFuncGaus->SetParNames("constant", "mean", "#sigma");
    hist->Fit("gaus", "RQM+", "", 30.0, 250.0);

    fitFuncLandau->GetParameters(&par[0]);
    fitFuncGaus->GetParameters(&par[3]);
    total->SetParameters(par);
    hist->Fit(total, "RQM+", "", 0.0, 250.0);

    delete fitFuncGaus;
    delete fitFuncLandau;
    return total;
  }
}

TF1 *Fits::Fit2Gaus(TH1D *hist) {
  if (hist->GetEntries() > 1000) {
    // if (hist->GetEntries() > 10000) ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    TF1 *fitFunc = new TF1("fitFunc", func::gausian2, min_value, max_value, 6);
    fitFunc->SetLineColor(2);
    par_max = std::isnan(hist->GetMaximum()) ? 0 : hist->GetMaximum();
    par_mean = std::isnan(hist->GetMean()) ? 0 : hist->GetMean();
    fitFunc->SetParameter(0, par_max);
    fitFunc->SetParameter(1, par_mean);
    fitFunc->SetParameter(2, 1);
    fitFunc->SetParNames("height", "mean", "FWHM");

    hist->Fit("fitFunc", "QM+", "", min_value, max_value);

    par_mean = std::isnan(fitFunc->GetParameter("mean")) ? 0 : fitFunc->GetParameter("mean");
    par_FWHM = std::isnan(fitFunc->GetParameter("FWHM")) ? 0 : fitFunc->GetParameter("FWHM");

    fitFunc->SetParameter(0, par_max);
    fitFunc->SetParameter(1, par_mean);
    fitFunc->SetParameter(2, par_FWHM);
    hist->Fit("fitFunc", "QM+", "", min_value, max_value);

    mean = fitFunc->GetParameter("mean");
    FWHM = fitFunc->GetParameter("FWHM");
    sigma = fitFunc->GetParameter("FWHM") / (2 * sqrt(2 * log(2)));  // 2.35482004503;
    return fitFunc;
  }
}

TF1 *Fits::FitLandau(TH1D *hist) {
  if (hist->GetEntries() > 1000) {
    // if (hist->GetEntries() > 10000) ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");

    TF1 *fitFunc = new TF1("fitFunc", "landau", -100.0, 100.0);
    fitFunc->SetLineColor(color);
    hist->Fit("fitFunc", "QM+", "", min_value, max_value);
    for (int i = 0; i < 10; i++) hist->Fit("fitFunc", "QM+", "", min_value, max_value);

    return fitFunc;
  }
}

TF1 *Fits::FitPoly_1D(TH1D *hist) {
  // if (hist->GetEntries() > 10000) ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  TF1 *fitFunc = new TF1("fitFunc", func::pol1, min_value, max_value);
  fitFunc->SetLineColor(9);
  fitFunc->SetParNames("intercept", "slope");

  hist->Fit("fitFunc", "QM+", "", min_value, max_value);

  fitFunc->SetParameter(0, fitFunc->GetParameter("intercept"));
  fitFunc->SetParameter(1, fitFunc->GetParameter("slope"));

  hist->Fit("fitFunc", "QM+", "", min_value, max_value);

  a = fitFunc->GetParameter("intercept");
  b = fitFunc->GetParameter("slope");
  return fitFunc;
}

TF1 *Fits::FitPoly_2D(TH1D *hist) {
  // if (hist->GetEntries() > 10000) ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  TF1 *fitFunc = new TF1("fitFunc", func::pol2, min_value, max_value);
  fitFunc->SetLineColor(30);
  fitFunc->SetParNames("a", "b", "c");

  hist->Fit("fitFunc", "QM+", "", min_value, max_value);

  fitFunc->SetParameter(0, fitFunc->GetParameter("a"));
  fitFunc->SetParameter(1, fitFunc->GetParameter("b"));
  fitFunc->SetParameter(2, fitFunc->GetParameter("c"));

  hist->Fit("fitFunc", "QM+", "", min_value, max_value);

  a = fitFunc->GetParameter("a");
  b = fitFunc->GetParameter("b");
  c = fitFunc->GetParameter("c");
  return fitFunc;
}

TF1 *Fits::FitPoly_3D(TH1D *hist) {
  // if (hist->GetEntries() > 10000) ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  TF1 *fitFunc = new TF1("fitFunc", func::pol3, min_value, max_value);
  fitFunc->SetLineColor(46);
  fitFunc->SetParNames("a", "b", "c", "d");

  hist->Fit("fitFunc", "QM+", "", min_value, max_value);

  fitFunc->SetParameter(0, fitFunc->GetParameter("a"));
  fitFunc->SetParameter(1, fitFunc->GetParameter("b"));
  fitFunc->SetParameter(2, fitFunc->GetParameter("c"));
  fitFunc->SetParameter(3, fitFunc->GetParameter("d"));

  hist->Fit("fitFunc", "QM+", "", min_value, max_value);

  a = fitFunc->GetParameter("a");
  b = fitFunc->GetParameter("b");
  c = fitFunc->GetParameter("c");
  d = fitFunc->GetParameter("d");
  return fitFunc;
}

TF1 *Fits::FitPoly_4D(TH1D *hist) {
  // if (hist->GetEntries() > 10000) ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  TF1 *fitFunc = new TF1("fitFunc", func::pol4, min_value, max_value);
  fitFunc->SetLineColor(42);
  fitFunc->SetParNames("a", "b", "c", "d", "e");

  hist->Fit("fitFunc", "QM+", "", min_value, max_value);

  fitFunc->SetParameter(0, fitFunc->GetParameter("a"));
  fitFunc->SetParameter(1, fitFunc->GetParameter("b"));
  fitFunc->SetParameter(2, fitFunc->GetParameter("c"));
  fitFunc->SetParameter(3, fitFunc->GetParameter("d"));
  fitFunc->SetParameter(4, fitFunc->GetParameter("e"));

  hist->Fit("fitFunc", "QM+", "", min_value, max_value);

  a = fitFunc->GetParameter("a");
  b = fitFunc->GetParameter("b");
  c = fitFunc->GetParameter("c");
  d = fitFunc->GetParameter("d");
  e = fitFunc->GetParameter("e");
  return fitFunc;
}

TF1 *Fits::FitPoly_fid(TH2D *hist) {
  // if (hist->GetEntries() > 10000) ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  TF1 *fitFunc = new TF1("fitFunc", "pol8", min_value, max_value);
  fitFunc->SetLineColor(42);
  fitFunc->SetParNames("a", "b", "c", "d", "e", "f", "g", "h");

  hist->Fit("fitFunc", "QM+", "", min_value, max_value);

  fitFunc->SetParameter(0, fitFunc->GetParameter("a"));
  fitFunc->SetParameter(1, fitFunc->GetParameter("b"));
  fitFunc->SetParameter(2, fitFunc->GetParameter("c"));
  fitFunc->SetParameter(3, fitFunc->GetParameter("d"));
  fitFunc->SetParameter(4, fitFunc->GetParameter("e"));
  fitFunc->SetParameter(5, fitFunc->GetParameter("f"));
  fitFunc->SetParameter(6, fitFunc->GetParameter("g"));
  fitFunc->SetParameter(7, fitFunc->GetParameter("h"));

  hist->Fit("fitFunc", "QM+", "", min_value, max_value);
  return fitFunc;
}

double Fits::fiducial_phi_lo(double theta_e, double theta_e_min, double k, double m, int c) {
  return -c * pow(sin((theta_e - theta_e_min) * 0.01745), k + m / theta_e + 1500. / (theta_e * theta_e));
}

double Fits::fiducial_phi_hi(double theta_e, double theta_e_min, double k, double m, int c) {
  return c * pow(sin((theta_e - theta_e_min) * 0.01745), k + m / theta_e + 1500. / (theta_e * theta_e));
}

TF1 *Fits::FitFiducial(TGraph *profile) {
  // ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  TF1 *fitFunc = new TF1("fitFunc", func::fiducial_phi, min_value, max_value, 10);

  fitFunc->SetLineColor(7);

  profile->Fit("fitFunc", "QM0+", "", min_value, max_value);

  fitFunc->SetParameter(0, fitFunc->GetParameter(0));
  fitFunc->SetParameter(1, fitFunc->GetParameter(1));
  fitFunc->SetParameter(2, fitFunc->GetParameter(2));
  fitFunc->SetParameter(3, fitFunc->GetParameter(3));
  fitFunc->SetParameter(4, fitFunc->GetParameter(4));
  fitFunc->SetParameter(5, fitFunc->GetParameter(5));
  fitFunc->SetParameter(6, fitFunc->GetParameter(6));
  fitFunc->SetParameter(7, fitFunc->GetParameter(7));
  fitFunc->SetParameter(8, fitFunc->GetParameter(8));
  fitFunc->SetParameter(9, fitFunc->GetParameter(9));
  fitFunc->SetParameter(10, fitFunc->GetParameter(10));

  profile->Fit("fitFunc", "QM+", "", min_value, max_value);

  return fitFunc;
}

TF1 *Fits::FitFiducial_hi(TH2D *hist2d) {
  // ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  a = b = c = d = 0.5;
  TF1 *fitFunc_hi = new TF1("fitFunc_hi",
                            "[0]*TMath::Power(TMath::Sin((x-[1])*"
                            "0.01745),[2] +[3] / x + 1500. / (x * x))",
                            max_value, max_value);

  fitFunc_hi->SetLineColor(8);
  fitFunc_hi->SetParNames("a", "b", "c", "d");

  hist2d->Fit("fitFunc_hi", "QM+", "", max_value, max_value);

  fitFunc_hi->SetParameter(0, fitFunc_hi->GetParameter("a"));
  fitFunc_hi->SetParameter(1, fitFunc_hi->GetParameter("b"));
  fitFunc_hi->SetParameter(2, fitFunc_hi->GetParameter("c"));
  fitFunc_hi->SetParameter(3, fitFunc_hi->GetParameter("d"));

  hist2d->Fit("fitFunc_hi", "QM+", "", max_value, max_value);

  a = fitFunc_hi->GetParameter("a");
  b = fitFunc_hi->GetParameter("b");
  c = fitFunc_hi->GetParameter("c");
  d = fitFunc_hi->GetParameter("d");
  return fitFunc_hi;
}

TF1 *Fits::FitFiducial(TH2D *hist2d) {
  // ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  a = b = c = d = 0.5;
  TF1 *fitFunc = new TF1("fitFunc",
                         "[0]*TMath::Power(TMath::Sin((x-[1])*"
                         "0.01745),[2] +[3] / x + 1500. / (x * x))",
                         max_value, max_value);

  fitFunc->SetLineColor(9);
  fitFunc->SetParNames("a", "b", "c", "d");

  hist2d->Fit("fitFunc", "QM+", "", max_value, max_value);

  fitFunc->SetParameter(0, fitFunc->GetParameter("a"));
  fitFunc->SetParameter(1, fitFunc->GetParameter("b"));
  fitFunc->SetParameter(2, fitFunc->GetParameter("c"));
  fitFunc->SetParameter(3, fitFunc->GetParameter("d"));

  hist2d->Fit("fitFunc", "QM+", "", max_value, max_value);

  a = fitFunc->GetParameter("a");
  b = fitFunc->GetParameter("b");
  c = fitFunc->GetParameter("c");
  d = fitFunc->GetParameter("d");
  return fitFunc;
}

TF1 *Fits::FitGenNormal(TH1D *hist) {
  TF1 *fitFunc = new TF1("genNormal", func::genNormal, min_value, max_value, 4);
  double min, max, val, min_m, max_m;
  if (hist->GetEntries() > 1000) {
    // if (hist->GetEntries() > 10000) ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");

    fitFunc->SetParLimits(1, 5.0, 200.0);

    fitFunc->SetParameter(0, 15.0);
    fitFunc->SetParameter(1, 10.0);
    fitFunc->SetParameter(2, (min_value + max_value) / 2.0);
    fitFunc->SetParameter(3, 3000.0);
    fitFunc->SetParNames("alpha", "beta", "mu", "weight");

    for (int i = 0; i < 10; i++) hist->Fit("genNormal", "QMR0+", "", min_value, max_value);
    hist->Fit("genNormal", "QMR+", "", min_value, max_value);

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
    left_edge_x = min_m;
    right_edge_x = max_m;
  }
  return fitFunc;
}

TF1 *Fits::FitBreitWigner(TH1D *hist) {
  // if (hist->GetEntries() > 10000) ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit");

  TF1 *fitbw = new TF1("bw", func::breit_wigner, min_value, max_value, 3);
  par_max = std::isnan(hist->GetMaximum()) ? 0 : hist->GetMaximum();
  par_mean = std::isnan(hist->GetMean()) ? 0 : hist->GetMean();
  par_RMS = std::isnan(hist->GetRMS()) ? 0 : hist->GetRMS();
  fitbw->SetParameter(0, par_mean);
  fitbw->SetParameter(1, par_RMS);
  fitbw->SetParameter(2, par_max);
  fitbw->SetParNames("Mean", "Width", "Const");

  for (int i = 0; i < 10; i++) {
    hist->Fit("bw", "QM0+", "", min_value, max_value);
  }
  hist->Fit("bw", "QM+", "", min_value, max_value);
  return fitbw;
}
