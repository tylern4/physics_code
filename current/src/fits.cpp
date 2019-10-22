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
  if (hist->GetEntries() > 10000) ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  if (hist->GetEntries() > 100) {
    TF1 *fitFunc = new TF1("fitFunc", func::gausian, -100.0, 100.0, 3);

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
  return NULL;
}

TF1 *Fits::FitLandauGaus(TH1D *hist) {
  if (hist->GetEntries() > 10000) ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  if (hist->GetEntries() > 1000) {
    double par[6];

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
  return NULL;
}

TF1 *Fits::Fit2Gaus(TH1D *hist) {
  if (hist->GetEntries() > 10000) ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  if (hist->GetEntries() > 1000) {
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
  return NULL;
}

TF1 *Fits::FitLandau(TH1D *hist) {
  if (hist->GetEntries() > 10000) ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  if (hist->GetEntries() > 1000) {
    TF1 *fitFunc = new TF1("fitFunc", "landau", -100.0, 100.0);
    fitFunc->SetLineColor(color);
    hist->Fit("fitFunc", "QM+", "", min_value, max_value);
    for (int i = 0; i < 10; i++) hist->Fit("fitFunc", "QM+", "", min_value, max_value);

    return fitFunc;
  }
  return NULL;
}

TF1 *Fits::FitPoly_1D(TH1D *hist) {
  if (hist->GetEntries() > 10000) ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");

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
  if (hist->GetEntries() > 10000) ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");

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
  if (hist->GetEntries() > 10000) ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");

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
  if (hist->GetEntries() > 10000) ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");

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

TF1 *Fits::FitPoly_fid(TGraph *hist) {
  TF1 *fitFunc = new TF1("fitFunc", func::pol4, min_value, max_value, 5);
  fitFunc->SetLineColor(46);

  hist->Fit("fitFunc", "QMR+", "", min_value, max_value);

  hist->Fit("fitFunc", "QMR+", "", min_value, max_value);
  return fitFunc;
}

double Fits::fiducial_phi_lo(double theta_e, double theta_e_min, double k, double m, int c) {
  return -c * pow(sin((theta_e - theta_e_min) * 0.01745), k + m / theta_e + 1500. / (theta_e * theta_e));
}

double Fits::fiducial_phi_hi(double theta_e, double theta_e_min, double k, double m, int c) {
  return c * pow(sin((theta_e - theta_e_min) * 0.01745), k + m / theta_e + 1500. / (theta_e * theta_e));
}

TF1 *Fits::FitFiducial(TGraph *profile, int sec) {
  TF1 *fitFunc = new TF1("fit_fid", func::fiducial, min_phi[sec], max_phi[sec], 6);
  fitFunc->SetParameters(1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
  fitFunc->SetParLimits(0, -10, 10);
  fitFunc->SetParLimits(1, -10, 10);
  fitFunc->SetParLimits(2, -10, 10);
  fitFunc->SetParLimits(3, -10, 10);
  fitFunc->SetParLimits(4, -10, 10);
  fitFunc->SetParLimits(5, -10, 10);

  fitFunc->SetLineColor(48);
  for (size_t i = 0; i < 10; i++) {
    fitFunc->SetParameter(0, fitFunc->GetParameter(0));
    fitFunc->SetParameter(1, fitFunc->GetParameter(1));
    fitFunc->SetParameter(2, fitFunc->GetParameter(2));
    fitFunc->SetParameter(3, fitFunc->GetParameter(3));
    fitFunc->SetParameter(4, fitFunc->GetParameter(4));
    fitFunc->SetParameter(5, fitFunc->GetParameter(5));

    profile->Fit("fit_fid", "QM0+", "", min_phi[sec], max_phi[sec]);
  }

  return fitFunc;
}

TF1 *Fits::FitFiducial_hi(TH2D *hist2d) {
  if (hist2d->GetEntries() > 10000) ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");

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
  if (hist2d->GetEntries() > 10000) ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");

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
  if (hist->GetEntries() > 10000) ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  TF1 *fitFunc = new TF1("genNormal", func::genNormal, min_value, max_value, 4);
  double min, max, val, min_m, max_m;
  if (hist->GetEntries() < 100) return NULL;

  fitFunc->SetParLimits(1, 2.0, 200.0);
  fitFunc->SetParameter(0, 15.0);
  fitFunc->SetParameter(1, 10.0);
  fitFunc->SetParameter(2, (min_value + max_value) / 2.0);
  fitFunc->SetParameter(3, 3000.0);
  fitFunc->SetParNames("alpha", "beta", "mu", "weight");

  for (int i = 0; i < 10; i++) hist->Fit("genNormal", "QMR0+", "", min_value, max_value);

  hist->Fit("genNormal", "QMR+", "", min_value, max_value);

  for (double m = min_value; m < max_value; m = m + 0.001) {
    // val = fitFunc->Derivative3(m);
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

  return fitFunc;
}

TF1 *Fits::FitBreitWigner(TH1D *hist) {
  if (hist->GetEntries() > 1000) ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");

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

TF1 *Fits::FitBreitWigner(std::shared_ptr<TH1D> &hists) { return Fits::FitBreitWigner(hists.get()); }

TF1 *Fits::FitMissMass(TH1D *hist) {
  if (hist->GetEntries() > 10000) ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  if (hist->GetEntries() < 1000) return nullptr;

  float fit_min = 0.5;
  float fit_max = 1.8;
  static const int max_par = 15;
  TF1 *total = new TF1("total", func::missMassfitFunction, fit_min, fit_max, max_par);
  TF1 *back_fit = new TF1("back_fit", func::missMassbackground, fit_min, fit_max, 6);
  TF1 *back_peak_fit = new TF1("back_peak_fit", func::missMasspeak, fit_min, fit_max, 3);
  back_peak_fit->SetParameters(1.0, 1.0, 1.0);

  TF1 *back_peak2_fit = new TF1("back_peak2_fit", func::missMasspeak, fit_min, fit_max, 3);
  back_peak2_fit->SetParameters(1.0, 1.0, 1.0);

  TF1 *peak_fit = new TF1("peak_fit", func::missMasspeak, fit_min, fit_max, 3);
  peak_fit->SetParameters(1.0, 1.0, 1.0);

  Double_t par[max_par];
  for (size_t i = 0; i < 10; i++) hist->Fit(peak_fit, "RNQM+", "", 0.9, 1.0);  // Peak of N at 0.939
  for (size_t i = 0; i < 50; i++) hist->Fit(back_fit, "RNQM+", "", 0.5, fit_max);
  for (size_t i = 0; i < 10; i++) hist->Fit(back_peak_fit, "RNQM+", "", 1.1, 1.3);
  for (size_t i = 0; i < 10; i++) hist->Fit(back_peak2_fit, "RNQM+", "", 1.3, 1.6);

  peak_fit->GetParameters(&par[0]);
  back_peak_fit->GetParameters(&par[3]);
  back_peak2_fit->GetParameters(&par[6]);
  back_fit->GetParameters(&par[9]);

  total->SetParameters(par);
  total->SetParName(0, "C_{N}");
  total->SetParName(1, "#Gamma_{N}");
  total->SetParName(2, "#mu_{N}");

  total->SetParName(3, "C_{#Delta^{0}}");
  total->SetParName(4, "#Gamma_{#Delta^{0}}");
  total->SetParName(5, "#mu_{#Delta^{0}}");

  total->SetParName(6, "C_{?}");
  total->SetParName(7, "#Gamma_{?}");
  total->SetParName(8, "#mu_{?}");

  total->SetParName(9, "C_{0}");
  total->SetParName(10, "C_{1}");
  total->SetParName(11, "C_{2}");
  total->SetParName(12, "C_{3}");
  total->SetParName(13, "C_{4}");
  total->SetParName(14, "C_{5}");

  for (size_t i = 0; i < 10; i++) hist->Fit(total, "RQNM+");

  hist->Fit(total, "RQM+");

  sigma = total->GetParameter("#Gamma_{N}") / (2 * sqrt(2 * log(2)));
  mean = total->GetParameter("#mu_{N}");

  return total;
}
