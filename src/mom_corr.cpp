#include "mom_corr.hpp"
#include <thread>
#include "color.hpp"

MomCorr::MomCorr() {
  std::cout << RED << "===== " << BLUE << "Using Momentum Corrections" << RED << " =====" << DEF << std::endl;
  if (getenv("MOM_CORR") != NULL) {
    _datadir = getenv("MOM_CORR");
    // Make a new thread to load each set of parameters
    std::thread theta_par(&MomCorr::read_theta_par, this);
    std::thread mom_par(&MomCorr::read_mom_par, this);
    std::thread mom_pip_par(&MomCorr::read_mom_pip_par, this);

    theta_par.join();
    mom_par.join();
    mom_pip_par.join();

  } else {
    std::cerr << "Set env MOM_CORR";
    exit(99);
  }
}

int MomCorr::GetSector(float phi) {
  // phi between -180 and 180
  // phi2 between 0 and 360
  float phi2 = phi;
  if (phi < 0) phi2 = phi2 + 360;
  int sect = 1 + ((int)phi2 + 30) / 60;
  if (sect == 7) sect = 1;
  return sect;
}

void MomCorr::read_theta_par() {
  /* Reading angle correction parameters */

  memset(&c0_theta[0][0], 0, ThetaC_n * NUM_SECTORS * sizeof(float));
  memset(&c1_theta[0][0], 0, ThetaC_n * NUM_SECTORS * sizeof(float));
  memset(&c2_theta[0][0], 0, ThetaC_n * NUM_SECTORS * sizeof(float));

  char file[100];
  for (int s = 1; s <= NUM_SECTORS; s++) {
    sprintf(file, "%s/angles_s%d.data", _datadir, s);
    FILE *fp = fopen(&file[0], "r");
    if (fp) {
      for (int j = 0; j < ThetaC_n; j++) {
        float x, y;
        float a1, a2, a3;
        char str[200];

        while (fgets(str, sizeof(str), fp)) {
          sscanf(str, "%f %f %f %f %f %f %f", &x, &a1, &y, &a2, &y, &a3, &y);
          int bin = (int)((x - ThetaC_min) / ThetaC_wid);
          c0_theta[bin][s - 1] = a1;
          c1_theta[bin][s - 1] = a2;
          c2_theta[bin][s - 1] = a3;
        }
      }
      fclose(fp);
    } else
      fprintf(stdout, "===>>> WARNING: Cannot read angle correction from file: %s\n", file);
  }
}

void MomCorr::read_mom_par() {
  /* Reading momentum correction parameters for electrons */

  memset(&c0_mom[0][0][0], 0, MomC_T_n * NUM_SECTORS * Npar * sizeof(float));
  memset(&c1_mom[0][0][0], 0, MomC_T_n * NUM_SECTORS * Npar * sizeof(float));

  char file[100];
  for (int s = 1; s <= NUM_SECTORS; s++) {
    for (int k = 0; k < 2; k++) {
      sprintf(file, "%s/momentum2_s%d_c%d.data", _datadir, s, k);
      FILE *fp = fopen(&file[0], "r");
      if (fp) {
        for (int j = 0; j < MomC_T_n; j++) {
          float x1, y;
          float a1, a2, a3, a4;
          char str[250];

          while (fgets(str, sizeof(str), fp)) {
            sscanf(str, "%f %f %f %f %f %f %f %f %f", &x1, &a1, &y, &a2, &y, &a3, &y, &a4, &y);

            int bin = (int)((x1 - MomC_T_min) / MomC_T_wid);
            if (k == 0) {
              c0_mom[bin][s - 1][0] = a1;
              c0_mom[bin][s - 1][1] = a2;
              c0_mom[bin][s - 1][2] = a3;
              c0_mom[bin][s - 1][3] = a4;
              // fprintf(stdout, "%6.2f     %12.8f  %12.8f  %12.8f  %12.8f \n", x1, c0_mom[bin][s - 1][0],
              //        c0_mom[bin][s - 1][1], c0_mom[bin][s - 1][2], c0_mom[bin][s - 1][3]);
            } else if (k == 1) {
              c1_mom[bin][s - 1][0] = a1;
              c1_mom[bin][s - 1][1] = a2;
              c1_mom[bin][s - 1][2] = a3;
              c1_mom[bin][s - 1][3] = a4;
              // fprintf(stdout, "%6.2f     %12.8f  %12.8f  %12.8f  %12.8f \n", x1, c1_mom[bin][s - 1][0],
              //        c1_mom[bin][s - 1][1], c1_mom[bin][s - 1][2], c1_mom[bin][s - 1][3]);
            }
          }
        }
        fclose(fp);
      } else
        fprintf(stdout,
                "===>>> WARNING: Cannot read electron momentum correction from "
                "file: %s\n",
                file);
    }
  }
}

void MomCorr::read_mom_pip_par() {
  memset(&d0_mom[0][0][0], 0, MomC_T_n * NUM_SECTORS * Npar * sizeof(float));
  memset(&d1_mom[0][0][0], 0, MomC_T_n * NUM_SECTORS * Npar * sizeof(float));

  char file[100];
  for (int s = 1; s <= NUM_SECTORS; s++) {
    for (int k = 0; k < 2; k++) {
      sprintf(file, "%s/momentum3_s%d_c%d.data", _datadir, s, k);
      FILE *fp = fopen(&file[0], "r");
      if (fp) {
        for (int j = 0; j < MomC_T_n; j++) {
          float x1, y;
          float a1, a2, a3;
          char str[250];

          while (fgets(str, sizeof(str), fp)) {
            sscanf(str, "%f %f %f %f %f %f %f", &x1, &a1, &y, &a2, &y, &a3, &y);

            int bin = (int)((x1 - MomC_T_min) / MomC_T_wid);
            if (k == 0) {
              d0_mom[bin][s - 1][0] = a1;
              d0_mom[bin][s - 1][1] = a2;
              d0_mom[bin][s - 1][2] = a3;
            } else if (k == 1) {
              d1_mom[bin][s - 1][0] = a1;
              d1_mom[bin][s - 1][1] = a2;
              d1_mom[bin][s - 1][2] = a3;
            }
          }
        }
        fclose(fp);
      } else
        std::cerr << "===>>> WARNING: Cannot read pi+ momentum correction from file:" << file << std::endl;
    }
  }
}

/* ================================================================= */
float MomCorr::theta_corr(float ThetaM, float PhiM, int Sector) {
  /* input: phi between 0 and 360 deg */

  float phis = PhiM - 60 * (Sector - 1);
  if (phis > 330.) phis = phis - 360;

  int bin = (int)((ThetaM - ThetaC_min) / ThetaC_wid);
  float ThetaBin = ThetaC_min + ThetaC_wid * (bin + 0.5);

  int bin2 = bin + 1;
  if (ThetaM < ThetaBin) bin2 = bin - 1;
  if (bin2 < 0) bin2 = 0;
  if (bin2 >= ThetaC_n) bin2 = ThetaC_n - 1;
  float ThetaBin2 = ThetaC_min + ThetaC_wid * (bin2 + 0.5);

  float ThetaC =
      ThetaM - (c0_theta[bin][Sector - 1] + c1_theta[bin][Sector - 1] * phis + c2_theta[bin][Sector - 1] * phis * phis);
  float ThetaC2 = ThetaM - (c0_theta[bin2][Sector - 1] + c1_theta[bin2][Sector - 1] * phis +
                            c2_theta[bin2][Sector - 1] * phis * phis);
  float fw = abs(ThetaM - ThetaBin) / ThetaC_wid;
  float fw2 = 1. - fw;
  float Theta = fw2 * ThetaC + fw * ThetaC2;

  return Theta;
}

float MomCorr::mom_corr(float MomM, float ThetaM, float PhiM, int Sector) {
  /* input: phi between 0 and 360 deg */
  float phis = PhiM - 60 * (Sector - 1);
  if (phis > 330.) phis = phis - 360;

  int bin = (int)((ThetaM - MomC_T_min) / MomC_T_wid);
  float ThetaBin = MomC_T_min + MomC_T_wid * (bin + 0.5);
  float fw = abs(ThetaM - ThetaBin) / MomC_T_wid;

  float A0 = 0.;
  float A1 = 0.;
  for (int j = 0; j < Npar; j++) {
    A0 = A0 + c0_mom[bin][Sector - 1][j] * pow(MomM, j);
    A1 = A1 + c1_mom[bin][Sector - 1][j] * pow(MomM, j);
  }
  float MomC = MomM - (A0 + A1 * phis);

  int bin2 = bin + 1;
  if (ThetaM < ThetaBin) bin2 = bin - 1;
  if (bin2 < 0) bin2 = 0;
  if (bin2 >= MomC_T_n) bin2 = MomC_T_n - 1;
  float ThetaBin2 = MomC_T_min + MomC_T_wid * (bin2 + 0.5);
  float fw2 = 1. - fw;

  A0 = A1 = 0.;
  for (int j = 0; j < Npar; j++) {
    A0 = A0 + c0_mom[bin2][Sector - 1][j] * pow(MomM, j);
    A1 = A1 + c1_mom[bin2][Sector - 1][j] * pow(MomM, j);
  }
  float MomC2 = MomM - (A0 + A1 * phis);

  float Mom = fw2 * MomC + fw * MomC2;

  return Mom;
}

float MomCorr::mom_corr_pip(float MomM, float ThetaM, float PhiM, int Sector) {
  /* input: phi between 0 and 360 deg */
  float phis = PhiM - 60 * (Sector - 1);
  if (phis > 330.) phis = phis - 360;

  int bin = (int)((ThetaM - MomC_T_min) / MomC_T_wid);
  float ThetaBin = MomC_T_min + MomC_T_wid * (bin + 0.5);
  float fw = abs(ThetaM - ThetaBin) / MomC_T_wid;

  float A0 = 0.;
  float A1 = 0.;
  for (int j = 0; j < Npar; j++) {
    A0 = A0 + d0_mom[bin][Sector - 1][j] * pow(MomM, j);
    A1 = A1 + d1_mom[bin][Sector - 1][j] * pow(MomM, j);
  }
  float MomC = MomM - (A0 + A1 * phis);

  int bin2 = bin + 1;
  if (ThetaM < ThetaBin) bin2 = bin - 1;
  if (bin2 < 0) bin2 = 0;
  if (bin2 >= MomC_T_n) bin2 = MomC_T_n - 1;
  float ThetaBin2 = MomC_T_min + MomC_T_wid * (bin2 + 0.5);
  float fw2 = 1. - fw;

  A0 = A1 = 0.;
  for (int j = 0; j < Npar; j++) {
    A0 = A0 + d0_mom[bin2][Sector - 1][j] * pow(MomM, j);
    A1 = A1 + d1_mom[bin2][Sector - 1][j] * pow(MomM, j);
  }
  float MomC2 = MomM - (A0 + A1 * phis);

  float Mom = fw2 * MomC + fw * MomC2;

  return Mom;
}

/* ================================================================= */
std::shared_ptr<LorentzVector> MomCorr::CorrectedVector(float px, float py, float pz, int particle_type) {
  auto Pin = std::make_unique<LorentzVector>(px, py, pz, mass_map[particle_type]);
  float theta = Pin->Theta() * RAD2DEG;
  float phi = Pin->Phi() * RAD2DEG;
  if (phi < 0) phi = phi + 360.;
  int sect = GetSector(phi);

  // if anything is bad then just return an uncorrected vector and let ROOT handle it
  // Only really effects MC events
  if (sect == 0 || sect > NUM_SECTORS || std::isnan(theta) || std::isnan(phi))
    return std::make_shared<LorentzVector>(px, py, pz, mass_map[particle_type]);

  float theta_c = theta_corr(theta, phi, sect);
  float mom_c = Pin->P();
  if (particle_type == ELECTRON)
    mom_c = mom_corr(Pin->P(), theta, phi, sect);
  else if (particle_type == PIM)
    mom_c = mom_corr(Pin->P(), theta, phi, sect);
  else if (particle_type == KM)
    mom_c = mom_corr(Pin->P(), theta, phi, sect);
  else if (particle_type == PIP)
    mom_c = mom_corr_pip(Pin->P(), theta, phi, sect);
  else if (particle_type == PROTON)
    mom_c = mom_corr_pip(Pin->P(), theta, phi, sect);
  else if (particle_type == KP)
    mom_c = mom_corr_pip(Pin->P(), theta, phi, sect);

  float _px = mom_c * sinf(theta_c * DEG2RAD) * cosf(Pin->Phi());
  float _py = mom_c * sinf(theta_c * DEG2RAD) * sinf(Pin->Phi());
  float _pz = mom_c * cosf(theta_c * DEG2RAD);

  return std::make_shared<LorentzVector>(_px, _py, _pz, mass_map[particle_type]);
}
