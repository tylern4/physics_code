/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef DATAHANDELER_H_GUARD
#define DATAHANDELER_H_GUARD
#include "main.h"
// Mashing together W vs Q2 and Delta T cuts into one file
// Saving the old files in a new folder to refer back to.
//
void dataHandeler(char *fin, char *RootFile_output, bool first_run) {
  TFile *RootOutputFile;
  // My Classes
  Histogram *hists = new Histogram();
  // From missing_mass.hpp :: missing_mass_calc()
  MissingMass *MM_neutron = new MissingMass();
  MM_neutron->Set_target_mass(MASS_P);
  MM_neutron->Set_target_PxPyPz(0);

  TCanvas *c1 = new TCanvas("c1", "c1", 100, 100);
  auto number_cols = 0;
  char rootFile[500];
  int num_of_events, total_events;
  bool cuts, electron_cuts;
  double e_E;

  std::vector<bool> pip_vec(MAX_PARTS, false);
  std::vector<bool> pim_vec(MAX_PARTS, false);
  std::vector<bool> proton_vec(MAX_PARTS, false);
  std::vector<bool> elec_vec(MAX_PARTS, false);

  // in main.h now
  // ofstream cut_outputs;
  // cut_outputs.open("outputFiles/cut_outputs.csv");
  // cut_outputs << "Cut,Mean,Sigma" << endl;

  int num_of_pis, num_of_proton;

  TVector3 e_mu_prime_3;
  TLorentzVector e_mu_prime;
  TLorentzVector e_mu(0.0, 0.0, sqrt(Square(E1D_E0) - Square(MASS_E)), E1D_E0);

  TVector3 Particle3(0.0, 0.0, 0.0);
  TLorentzVector Particle4(0.0, 0.0, 0.0, 0.0);
  double theta, phi;
  int sector;
  // End declrare variables

  // Open outputfile
  RootOutputFile = new TFile(RootFile_output, "RECREATE");

  // Load chain from branch h10
  TChain chain("h10");
  TString fins = fin;
  TFileCollection fc("fileList", "", fins.Data());
  chain.AddFileInfoList((TCollection *)fc.GetList());
  cout << blue << "Analyzing file " << green << fin << def << bgdef << endl;

  getBranches(&chain);
  if (!first_run)
    getMorebranchs(&chain);

  num_of_events = (int)chain.GetEntries();

  //#pragma omp parallel
  //{
  for (int current_event = 0; current_event < num_of_events; current_event++) {
    // update loadbar and get current event
    loadbar(current_event + 1, num_of_events);
    chain.GetEntry(current_event);

    // reset electron cut bool
    electron_cuts = true;
    // electron cuts
    electron_cuts &= (ec[0] > 0); // ``` ``` ``` ec
    electron_cuts &= ((int)id[0] == ELECTRON ||
                      (int)id[0] == 0); // First particle is electron
    electron_cuts &=
        ((int)gpart > 0); // Number of good particles is greater than 0
    electron_cuts &= ((int)stat[0] > 0); // First Particle hit stat
    electron_cuts &= ((int)q[0] == -1);  // First particle is negative Q
    electron_cuts &= ((int)sc[0] > 0);   // First Particle hit sc
    electron_cuts &= ((int)dc[0] > 0);   // ``` ``` ``` dc
    electron_cuts &= ((int)dc_stat[dc[0] - 1] > 0);

    if (electron_cuts)
      hists->EC_fill(etot[ec[0] - 1], p[0]);

    electron_cuts &= (p[0] > MIN_P_CUT); // Minimum Momentum cut

    if (electron_cuts && cc[0] > 0) {
      int cc_sector = cc_sect[cc[0] - 1];
      int cc_segment = (cc_segm[0] % 1000) / 10;
      int cc_pmt = cc_segm[0] / 1000 - 1;
      int cc_nphe = nphe[cc[0] - 1];
      double theta_cc =
          TMath::ACos(TMath::Abs(p[0] * cz[0]) / TMath::Abs(p[0]));

      theta_cc = theta_cc / D2R;
      hists->CC_fill(cc_sector, cc_segment, cc_pmt, cc_nphe, theta_cc);

      hists->Fill_Beam_Position((double)dc_vx[dc[0] - 1],
                                (double)dc_vy[dc[0] - 1]);
    }

    if (electron_cuts) {
      if (first_run) {
        is_electron = &elec_vec;
        is_electron->at(0) = true;
        for (int part_num = 1; part_num < gpart; part_num++) {
          is_pip = &pip_vec;
          is_pim = &pim_vec;
          is_proton = &proton_vec;
          is_pip->at(part_num) = (id[part_num] == PIP);
          is_proton->at(part_num) = (id[part_num] == PROTON);
          is_pim->at(part_num) = (id[part_num] == PIM);
        }
      }

      // Setup scattered electron 4 vector
      e_mu_prime_3.SetXYZ(p[0] * cx[0], p[0] * cy[0], p[0] * cz[0]);
      e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);
      // Set the vertex time (time of electron hit)
      Delta_T *delta_t = new Delta_T();
      delta_t->delta_t_cut(hists, first_run);

      theta = theta_calc(cz[0]);
      // phi = center_phi_calc(cx[0],cy[0]);
      phi = phi_calc(cx[0], cy[0]);
      sector = get_sector(phi);

      // Fill_fid(theta,phi,get_sector(phi_calc(cx[0],cy[0])));

      hists->Fill_fid(theta, phi, sector);

      if (first_run) {
        W = W_calc(e_mu, e_mu_prime);
        Q2 = Q2_calc(e_mu, e_mu_prime);
        e_E = e_mu_prime.E();
      }

      hists->WvsQ2_Fill(e_E, W, Q2, xb_calc(Q2, e_E));
      num_of_proton = num_of_pis = 0;

      //#pragma omp parallel for
      for (int part_num = 1; part_num < gpart; part_num++) {
        if (p[part_num] == 0)
          continue;

        // if(is_proton->at(part_num) == is_pip->at(part_num)) continue;

        hists->Fill_Mass(m[part_num]);
        Particle3.SetXYZ(p[part_num] * cx[part_num], p[part_num] * cy[part_num],
                         p[part_num] * cz[part_num]);
        Particle4.SetVectM(Particle3, Get_Mass(id[part_num]));

        hists->MomVsBeta_Fill(Particle4.E(), p[part_num], b[part_num]);
        if (q[part_num] == 1) {
          hists->MomVsBeta_Fill_pos(p[part_num], b[part_num]);
          if (is_proton->at(part_num) && (id[part_num] == PROTON)) {
            num_of_proton++;
            hists->Fill_proton_WQ2(W, Q2);
            hists->Fill_proton_ID_P(p[part_num], b[part_num]);
          }
          if (is_pip->at(part_num) && (id[part_num] == PIP)) {
            num_of_pis++;
            hists->Fill_pion_WQ2(W, Q2);
            hists->Fill_Pi_ID_P(p[part_num], b[part_num]);
            TLorentzVector gamma_mu = (e_mu - e_mu_prime);
            if (first_run) {
              MM_neutron->Set_PxPyPz(p[part_num] * cx[part_num],
                                     p[part_num] * cy[part_num],
                                     p[part_num] * cz[part_num]);
              MM = MM_neutron->missing_mass(gamma_mu);
            }
          }

          if ((is_pip->at(part_num) && (id[part_num] == PIP)) ||
              (is_proton->at(part_num) && (id[part_num] == PROTON))) {
            hists->Fill_proton_Pi_ID_P(p[part_num], b[part_num]);
          }
        } else if (q[part_num] == -1) {
          hists->MomVsBeta_Fill_neg(p[part_num], b[part_num]);
        }
      }

      if (num_of_pis == 1)
        hists->Fill_single_pi_WQ2(W, Q2);

      if (num_of_pis == 1 && num_of_proton == 0) {
        hists->Fill_Missing_Mass(MM);
        hists->Fill_Missing_Mass_square(Square(MM));
      }
      if (num_of_proton == 1)
        hists->Fill_single_proton_WQ2(W, Q2);
    }
  }
  //}
  std::cout << green << "Fitting" << def << std::endl;
  // Start of cuts
  Fits MM_neutron_cut;
  double fit_range_min = 0.88;
  double fit_range_max = 1.0;
  MM_neutron_cut.FitGaus(hists->Missing_Mass, fit_range_min, fit_range_max);

  Header *MM_header = new Header("../src/missing_mass_gaussians.hpp", "MM");
  MM_header->WriteGaussian("mm", 1, MM_neutron_cut.mean, MM_neutron_cut.sigma);

  Fits MissingMassSquare_cut;
  fit_range_min = 0.5;
  fit_range_max = 1.1;
  MissingMassSquare_cut.FitGaus(hists->Missing_Mass_square, fit_range_min,
                                fit_range_max);
  MM_header->WriteGaussian("mm_square", 1, MissingMassSquare_cut.mean,
                           MissingMassSquare_cut.sigma);
  delete MM_header;

  //
  // end stuff
  chain.Reset(); // delete Tree object

  RootOutputFile->cd();
  hists->EC_Write();
  hists->Beam_Position_Write();

  TDirectory *WvsQ2_folder = RootOutputFile->mkdir("W vs Q2");
  WvsQ2_folder->cd();
  hists->WvsQ2_Write();

  TDirectory *MomVsBeta_folder = RootOutputFile->mkdir("Momentum vs beta");
  MomVsBeta_folder->cd();
  hists->MomVsBeta_Write();

  // Missing Mass Write
  TDirectory *MissMass = RootOutputFile->mkdir("Missing_Mass");
  MissMass->cd();
  hists->Write_Missing_Mass();

  // Delta T Write
  TDirectory *DeltaT = RootOutputFile->mkdir("Delta_T");
  DeltaT->cd();
  hists->delta_t_Write();

  TDirectory *DeltaT_slices = RootOutputFile->mkdir("Delta_T_slices");
  DeltaT_slices->cd();
  hists->delta_t_slices_Write();

  TDirectory *DeltaT_sec_pad = RootOutputFile->mkdir("Delta_T_sec_pad");
  DeltaT_sec_pad->cd();
  hists->delta_t_sec_pad_Write();

  TDirectory *Delta_T_canvases = RootOutputFile->mkdir("Delta_T_canvases");
  Delta_T_canvases->cd();
  hists->delta_T_canvas();

  TDirectory *Theta_CC_hists = RootOutputFile->mkdir("Theta_CC_hists");
  Theta_CC_hists->cd();
  hists->Theta_CC_Write();

  TDirectory *CC_hists = RootOutputFile->mkdir("CC_hists");
  CC_hists->cd();
  hists->CC_Write();

  TDirectory *CC_canvases = RootOutputFile->mkdir("CC_canvases");
  CC_canvases->cd();
  hists->CC_canvas();

  TDirectory *Fid_cuts = RootOutputFile->mkdir("Fid_cuts");
  Fid_cuts->cd();
  hists->Fid_Write();

  delete hists;
  RootOutputFile->Write();
  RootOutputFile->Close();
  cut_outputs.close();
}
#endif
