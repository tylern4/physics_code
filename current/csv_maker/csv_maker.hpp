/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef CSV_MAKER_H_GUARD
#define CSV_MAKER_H_GUARD
#include "main.h"
// Mashing together W vs Q2 and Delta T cuts into one file
// Saving the old files in a new folder to refer back to.
//
void make_electron_csv(char *fin, char *csv_name_output) {
  int num_of_events;
  bool electron_cuts;
  double e_E = 0.0;
  // in main.h now
  // ofstream cut_outputs;
  csv_output.open(csv_name_output);
  csv_output << "id_0, gpart, ec_0, stat_0, q_0, sc_0, dc_0, dc_stat, p, "
                "cx, cy, cz, cc_0, cc_sector, cc_segment, cc_pmt, cc_nphe, "
                "theta_cc, theta_fid, phi_fid, sector_fid, W, Q2, Energy"
             << endl;

  // Load chain from branch h10
  TChain chain("h10");

  chain.Add(fin);

  getBranches(&chain);

  num_of_events = (int)chain.GetEntries();

  //#pragma omp parallel for
  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain.GetEntry(current_event);
    electron_cuts = true;
    // electron cuts
    electron_cuts &= ((int)ec[0] > 0); // ``` ``` ``` ec
    electron_cuts &= ((int)id[0] == ELECTRON ||
                      (int)id[0] == 0);  // First particle is electron`
    electron_cuts &= ((int)gpart > 0);   // Number of good particles is gt 0
    electron_cuts &= ((int)stat[0] > 0); // First Particle hit stat
    electron_cuts &= ((int)q[0] == -1);  // First particle is negative Q
    electron_cuts &= ((int)sc[0] > 0);   // First Particle hit sc
    electron_cuts &= ((int)dc[0] > 0);   // ``` ``` ``` dc
    electron_cuts &= ((int)dc_stat[dc[0] - 1] > 0);

    csv_output << (int)id[0] << ",";
    csv_output << (int)gpart << ",";
    csv_output << (int)ec[0] << ",";
    csv_output << (int)stat[0] << ",";
    csv_output << (int)q[0] << ",";
    csv_output << (int)sc[0] << ",";
    csv_output << (int)dc[0] << ",";
    csv_output << (int)dc_stat[dc[0] - 1] << ",";
    csv_output << (double)p[0] << ",";
    csv_output << (double)cx[0] << ",";
    csv_output << (double)cy[0] << ",";
    csv_output << (double)cz[0] << ",";

    if (electron_cuts)
      electron_cuts &= (p[0] > MIN_P_CUT); // Minimum Momentum cut

    int cc_sector = 0;
    int cc_segment = 0;
    int cc_pmt = 0;
    int cc_nphe = 0;
    double theta_cc = 0.0;
    double theta = 0.0;
    double phi = 0.0;
    int sector = 0;

    if (electron_cuts && (int)cc[0] > 0) {
      cc_sector = cc_sect[(int)cc[0] - 1];
      cc_segment = ((int)cc_segm[0] % 1000) / 10;
      cc_pmt = (int)cc_segm[0] / 1000 - 1;
      cc_nphe = (int)nphe[(int)cc[0] - 1];
      theta_cc = TMath::ACos(TMath::Abs(p[0] * cz[0]) / TMath::Abs(p[0]));
      theta_cc = theta_cc / D2R;
    }
    csv_output << (int)cc[0] << ",";
    csv_output << cc_sector << ",";
    csv_output << cc_segment << ",";
    csv_output << cc_pmt << ",";
    csv_output << cc_nphe << ",";
    csv_output << theta_cc << ",";

    // Setup scattered electron 4 vector
    TVector3 e_mu_prime_3;
    TLorentzVector e_mu_prime;
    TLorentzVector e_mu(0.0, 0.0, sqrt(Square(E1D_E0) - Square(MASS_E)),
                        E1D_E0);

    e_mu_prime_3.SetXYZ(p[0] * cx[0], p[0] * cy[0], p[0] * cz[0]);
    e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);

    theta = theta_calc(cz[0]);
    phi = phi_calc(cx[0], cy[0]);
    sector = get_sector(phi);

    csv_output << theta << ",";
    csv_output << phi << ",";
    csv_output << sector << ",";

    W = W_calc(e_mu, e_mu_prime);
    Q2 = Q2_calc(e_mu, e_mu_prime);
    e_E = e_mu_prime.E();
    csv_output << W << ",";
    csv_output << Q2 << ",";
    csv_output << e_E;
    csv_output << endl;
  }
  //
  // end stuff
  chain.Reset();
  csv_output.close();
}
#endif
