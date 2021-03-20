#include "TCandle.h"
#include "TF1.h"
#include "TFile.h"
#include "TH2.h"
#include "TROOT.h"
#include "TStyle.h"

double thmax(double *x, double *par) {
  double s, y;
  s = 222.947 - 505.444 * x[0] + 534.834 * x[0] * x[0] - 223.395 * x[0] * x[0] * x[0];
  return s;
}

double th_s1_up1(double *x, double *par) {
  double s;

  s = (304.23 * (x[0] + 0.15) * (x[0] + 0.15) * (x[0] + 0.15) - 255.798 * (x[0] + 0.15) * (x[0] + 0.15) +
       497.462 * (x[0] + 0.15) + 38.0385) *
          exp(-1.85 * (x[0] + 0.15)) +
      5.5;

  return s;
}

double th_s1_down1(double *x, double *par) {
  double s;

  s = (304.23 * (x[0] + 0.15) * (x[0] + 0.15) * (x[0] + 0.15) - 255.798 * (x[0] + 0.15) * (x[0] + 0.15) +
       497.462 * (x[0] + 0.15) + 38.0385) *
          exp(-1.85 * (x[0] + 0.15)) -
      3;

  return s;
}

double th_s1_down2(double *x, double *par) {
  double s;

  s = (pow((x[0] - 0.103718), (0.0703664)) * 252.822 - 133.024) * exp(-0.5 * x[0]) + 0.1;

  return s;
}

double th_s1_up2(double *x, double *par) {
  double s;

  s = (pow((x[0] - 0.0575818), (0.075643)) * 238.248 - 115.039) * exp(-0.5 * x[0]) - 0.1;

  return s;
}

double th_s1_up3(double *x, double *par) {
  double s, y;

  y = x[0] + 0.25;

  s = (304.23 * y * y * y - 255.798 * y * y + 497.462 * y + 38.0385) * exp(-1.6 * y) - 104;

  return s;
}
double th_s1_down3(double *x, double *par) {
  double s, y;

  y = x[0] - 0.12;

  s = (304.23 * y * y * y - 255.798 * y * y + 497.462 * y + 38.0385) * exp(-1.6 * y) - 103;

  return s;
}
///-----------------------

double th_s2_up1(double *x, double *par) {
  double s;

  s = pow((x[0] - 0.415068), (0.226449)) * 48.7564 + 2.79478 - 1.;

  return s;
}

double th_s2_down1(double *x, double *par) {
  double s;

  s = pow((x[0] - 0.449975), (0.315164)) * 36.608 + 9.74262 - 1.;

  return s;
}

double th_s2_up2(double *x, double *par) {
  double s;

  s = (387.289 * x[0] * x[0] * x[0] - 758.466 * x[0] * x[0] + 842.881 * x[0] - 299.953 + 15.) * exp(-2 * x[0]);

  return s;
}

double th_s2_down2(double *x, double *par) {
  double s, y;
  y = x[0] + 0.03;

  s = (387.289 * y * y * y - 758.466 * y * y + 842.881 * y - 299.953 - 15.) * exp(-2 * y) - 1.5;
  return s;
}
//---------------------------------

double th_s3_up1(double *x, double *par) {
  double s;

  s = (10000 * x[0] * x[0] * x[0] - 3607.41 * x[0] * x[0] + 1725.72 * x[0] - 10.6776) * exp(-4.7 * x[0]);

  return s;
}

double th_s3_down1(double *x, double *par) {
  double s;

  s = (10000 * x[0] * x[0] * x[0] - 4505.62 * x[0] * x[0] + 2056.24 * x[0] - 77.4077 + 5.) * exp(-4.7 * x[0]);

  return s;
}

double th_s3_up2(double *x, double *par) {
  double s;

  s = pow((x[0] - 0.416536), (0.108376)) * 67.4593 - 21.4374;

  return s;
}

double th_s3_down2(double *x, double *par) {
  double s;

  s = pow((x[0] - 0.454898), (0.289291)) * 35.7267 + 6.65908 + 1.5;

  return s;
}
//-----------------------------------

double th_s4_up1(double *x, double *par) {
  double s, y;

  y = x[0] + 0.165;

  s = (304.23 * y * y * y - 255.798 * y * y + 497.462 * y + 38.0385) * exp(-1.85 * y) + 5.;

  return s;
}

double th_s4_down1(double *x, double *par) {
  double s;

  s = (304.23 * (x[0] + 0.18) * (x[0] + 0.18) * (x[0] + 0.18) - 255.798 * (x[0] + 0.18) * (x[0] + 0.18) +
       497.462 * (x[0] + 0.18) + 38.0385) *
          exp(-1.85 * (x[0] + 0.18)) -
      1.;

  return s;
}

double th_s4_up2(double *x, double *par) {
  double s;

  s = (1600 * (x[0] + 0.03) * (x[0] + 0.03) * (x[0] + 0.03) - 1068.36 * (x[0] + 0.03) * (x[0] + 0.03) +
       775.016 * (x[0] + 0.03) - 1.13034) *
      exp(-2.75 * (x[0] + 0.03));

  return s;
}

double th_s4_down2(double *x, double *par) {
  double s;

  s = (pow((x[0] - 0.103718), (0.0703664)) * 252.822 - 133.024) * exp(-0.45 * x[0]) - 7.;

  return s;
}
//---------------------------------
double th_s5_up1(double *x, double *par) {
  double s;

  s = (525.498 * (x[0] + 0.03) * (x[0] + 0.03) * (x[0] + 0.03) - 1284.98 * (x[0] + 0.03) * (x[0] + 0.03) +
       1460.67 * (x[0] + 0.03) - 499.999) *
      exp(-1.94 * (x[0] + 0.03));

  return s;
}

double th_s5_down1(double *x, double *par) {
  double s;

  s = (525.498 * (x[0] - 0.02) * (x[0] - 0.02) * (x[0] - 0.02) - 1284.98 * (x[0] - 0.02) * (x[0] - 0.02) +
       1460.67 * (x[0] - 0.02) - 499.999) *
          exp(-1.94 * (x[0] - 0.02)) -
      4.7;

  return s;
}

double th_s5_up2(double *x, double *par) {
  double s;

  s = (304.23 * (x[0]) * (x[0]) * (x[0]) - 255.798 * (x[0]) * (x[0]) + 497.462 * (x[0]) + 38.0385) *
      exp(-1.85 * (x[0]));

  return s;
}

double th_s5_down2(double *x, double *par) {
  double s;

  s = (304.23 * (x[0] + 0.03) * (x[0] + 0.03) * (x[0] + 0.03) - 255.798 * (x[0] + 0.03) * (x[0] + 0.03) +
       497.462 * (x[0] + 0.03) + 38.0385) *
          exp(-1.85 * (x[0] + 0.03)) -
      11.;

  return s;
}

double th_s5_up3(double *x, double *par) {
  double s;
  s = pow((x[0] - 0.304992), (0.0758186)) * 91.5643 - 48.2057 - 1.;
  return s;
}

double th_s5_down3(double *x, double *par) {
  double s;

  s = pow((x[0] - 0.36848), (0.0864219)) * 70.4769 - 34.9998 + 1.5;

  return s;
}
//-----------------------------------------------

double th_s6_up1(double *x, double *par) {
  double s;
  s = pow((x[0] - 0.05 - 0.0942469), (0.0582707)) * 114.358 - 50 - 0.5;
  return s;
}

double th_s6_down1(double *x, double *par) {
  double s;
  s = pow((x[0] - 0.05 - 0.126994), (0.0706829)) * 110.073 - 50 + 2.;
  return s;
}

double th_s6_down2(double *x, double *par) {
  double s;
  s = pow((x[0] - 0.454098), (0.0912936)) * 58.2946 - 20.4843 + 1.5;
  return s;
}

double th_s6_up2(double *x, double *par) {
  double s;
  s = pow((x[0] - 0.416536), (0.108376)) * 67.4593 - 21.4374 - 1.;
  return s;
}

double fit1(double *x, double *par) {
  double s, y;

  y = x[0] + par[0];
  s = (par[1] * y * y * y + par[2] * y * y + par[3] * y + par[4]) * exp(par[5] * y) + par[6];
  return s;
}

int th_vs_p_pip_all(std::string dumb) {
  TFile *f = new TFile(dumb.c_str(), "READ");

  TH2D *hist = (TH2D *)f->Get("Fid_cuts/pip_theta_p_nocut_1");
  TF1 *_fit = new TF1("fit1", fit1, 0, 4, 7);

  _fit->SetParameter(1, 304.23);
  _fit->SetParameter(2, -255.798);
  _fit->SetParameter(3, 497.462);
  _fit->SetParameter(4, 38.0385);
  _fit->SetParameter(5, -1.85);
  _fit->SetParameter(6, 11);
  hist->Fit(_fit);
  hist->Draw();

  // delete f;
  return 0;
};

#if not defined(__CLING__)
int main(int argc, char const *argv[]) {
  if (argc < 2) {
    std::cerr << "Not enough arguments" << std::endl;
    std::cerr << "To Use:\t" << argv[0] << " norad.root rad.root" << std::endl;
    exit(1);
  }

  return th_vs_p_pip_all(argv[1]);
}
#endif
//-----------------------------
