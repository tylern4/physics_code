/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#include "datahandeler.hpp"

DataHandeler::DataHandeler(std::vector<std::string> fin, std::string RootFile_output) {
  input_files = fin;
  pip_vec = new std::vector<bool>(MAX_PARTS, false);
  pim_vec = new std::vector<bool>(MAX_PARTS, false);
  proton_vec = new std::vector<bool>(MAX_PARTS, false);
  elec_vec = new std::vector<bool>(MAX_PARTS, false);
  c1 = new TCanvas("c1", "c1", 100, 100);

  // in main.h now
  // ofstream cut_outputs;
  // cut_outputs.open("outputFiles/cut_outputs.csv");
  // cut_outputs << "Cut,Mean,Sigma" << endl;

  e_mu = new TLorentzVector(0.0, 0.0, sqrt(E1D_E0 * E1D_E0 - MASS_E * MASS_E), E1D_E0);
  // End declrare variables

  // Open outputfile
  RootOutputFile = new TFile(RootFile_output.c_str(), "RECREATE");
  hists = new Histogram();
  MM_neutron = new MissingMass(MASS_P, 0.0);
}

DataHandeler::~DataHandeler() {
  delete MM_neutron;
  std::cout << GREEN << "\nFitting" << DEF << std::endl;
  // Start of cuts
  Fits *MM_neutron_cut = new Fits();
  MM_neutron_cut->Set_min(0.88);
  MM_neutron_cut->Set_max(1.0);
  // MM_neutron_cut->FitBreitWigner(hists->Missing_Mass);
  // MM_neutron_cut->FitGaus(hists->Missing_Mass);
  // MM_neutron_cut->Fit2Gaus(hists->Missing_Mass);
  MM_neutron_cut->FitLandau(hists->Missing_Mass);

  Header *MM_header = new Header("../src/missing_mass_gaussians.hpp", "MM");
  MM_header->WriteGaussian("mm", 1, MM_neutron_cut->Get_mean(), MM_neutron_cut->Get_sigma());

  Fits *MissingMassSquare_cut = new Fits();
  MissingMassSquare_cut->Set_min(0.5);
  MissingMassSquare_cut->Set_max(1.1);
  // MissingMassSquare_cut->FitBreitWigner(hists->Missing_Mass_square);
  // MissingMassSquare_cut->FitGaus(hists->Missing_Mass_square);
  // MissingMassSquare_cut->Fit2Gaus(hists->Missing_Mass_square);
  MissingMassSquare_cut->FitLandau(hists->Missing_Mass_square);
  MM_header->WriteGaussian("mm_square", 1, MissingMassSquare_cut->Get_mean(), MissingMassSquare_cut->Get_sigma());
  delete MM_header;
  delete MM_neutron_cut;

  //
  // end stuff

  RootOutputFile->cd();
  TDirectory *EC_folder = RootOutputFile->mkdir("EC_hists");
  EC_folder->cd();
  hists->EC_Write();

  TDirectory *EC_slices = RootOutputFile->mkdir("EC_slices");
  EC_slices->cd();
  hists->EC_slices_Write();

  TDirectory *Beam_Folder = RootOutputFile->mkdir("Beam Position");
  Beam_Folder->cd();
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

  TDirectory *Fid_canvas = RootOutputFile->mkdir("Fid_canvas");
  Fid_canvas->cd();
  hists->fid_canvas();

  delete hists;
  RootOutputFile->Write();
  RootOutputFile->Close();
  cut_outputs.close();
}

void DataHandeler::loadbar(long x, long n) {
  int w = 50;
  if ((x != n) && (x % (n / 100 + 1) != 0)) return;

  double ratio = x / (double)n;
  int c = ratio * w;

  std::cout << BLUE << " [";
  for (int x = 0; x < c; x++) std::cout << GREEN << "=" << DEF;
  std::cout << GREEN << ">" << DEF;
  for (int x = c; x < w; x++) std::cout << " ";
  std::cout << BLUE << (int)(ratio * 100) << "%]\r" << DEF << std::flush;
}

void DataHandeler::run() {
  int size = input_files.size();
  int i = 0;

#pragma omp parallel for private(i)
  for (i = 0; i < size; i++) {
    // loadbar(i + 1, size);
    file_handeler(input_files.at(i));
  }
}

void DataHandeler::file_handeler(std::string fin) {
  // Load chain from branch h10
  bool cuts, electron_cuts;
  int num_of_events;
  int total_events;
  int num_of_pis;
  int num_of_proton;
  double e_E;
  double theta;
  double phi;
  int sector;
  bool first_run = true;
  TChain *chain = new TChain("h10");
  chain->Add(fin.c_str());

  getBranches(chain);
  if (!first_run) getMorebranchs(chain);
  num_of_events = (int)chain->GetEntries();
  std::cout << "Hi I'm parallel " << num_of_events << std::endl;
  /*
    int current_event = 0;
    for (current_event = 0; current_event < num_of_events; current_event++) {
      chain->GetEntry(current_event);

      // reset electron cut bool
      electron_cuts = true;
      // electron cuts
      electron_cuts &= (ec[0] > 0);                                  // ``` ``` ``` ec
      electron_cuts &= ((int)id[0] == ELECTRON || (int)id[0] == 0);  // First particle is electron
      electron_cuts &= ((int)gpart > 0);                             // Number of good particles is greater than 0
      electron_cuts &= ((int)stat[0] > 0);                           // First Particle hit stat
      electron_cuts &= ((int)q[0] == -1);                            // First particle is negative Q
      electron_cuts &= ((int)sc[0] > 0);                             // First Particle hit sc
      electron_cuts &= ((int)dc[0] > 0);                             // ``` ``` ``` dc
      electron_cuts &= ((int)dc_stat[dc[0] - 1] > 0);
      electron_cuts &= ((int)cc[0] > 0);
      if (electron_cuts) hists->EC_fill(etot[ec[0] - 1], p[0]);
      electron_cuts &= (p[0] > MIN_P_CUT);  // Minimum Momentum cut????

      if (electron_cuts) {
        int cc_sector = cc_sect[cc[0] - 1];
        int cc_segment = (cc_segm[0] % 1000) / 10;
        int cc_pmt = cc_segm[0] / 1000 - 1;
        int cc_nphe = nphe[cc[0] - 1];
        double theta_cc = TMath::ACos(TMath::Abs(p[0] * cz[0]) / TMath::Abs(p[0]));

        theta_cc = theta_cc / D2R;
        hists->CC_fill(cc_sector, cc_segment, cc_pmt, cc_nphe, theta_cc);

        hists->Fill_Beam_Position((double)dc_vx[dc[0] - 1], (double)dc_vy[dc[0] - 1], (double)dc_vz[dc[0] - 1]);
        if (first_run) {
          is_electron = elec_vec;
          is_electron->at(0) = true;
          for (int part_num = 1; part_num < gpart; part_num++) {
            is_pip = pip_vec;
            is_pim = pim_vec;
            is_proton = proton_vec;
            is_pip->at(part_num) = (id[part_num] == PIP);
            is_proton->at(part_num) = (id[part_num] == PROTON);
            is_pim->at(part_num) = (id[part_num] == PIM);
          }
        }

        // Setup scattered electron 4 vector
        TLorentzVector e_mu_prime = physics::fourVec(p[0], cx[0], cy[0], cz[0], MASS_E);

        // Set the vertex time (time of electron hit)
        Delta_T *dt = new Delta_T(sc_t[sc[0] - 1], sc_r[sc[0] - 1]);
        dt->delta_t_hists(hists);
        // std::vector<double> dt_proton = dt->delta_t_array(MASS_P, gpart);
        // std::vector<double> dt_pi = dt->delta_t_array(MASS_PIP, gpart);
        delete dt;

        if (electron_cuts) {
          theta = physics::theta_calc(cz[0]);
          phi = physics::phi_calc(cx[0], cy[0]);
          sector = physics::get_sector(phi);

          hists->Fill_electron_fid(theta, phi, sector);
        }
        if (first_run) {
          W = physics::W_calc(*e_mu, e_mu_prime);
          Q2 = physics::Q2_calc(*e_mu, e_mu_prime);
          e_E = e_mu_prime.E();
        }

        hists->WvsQ2_Fill(e_E, W, Q2, physics::xb_calc(Q2, e_E));
        num_of_proton = num_of_pis = 0;
        for (int part_num = 1; part_num < gpart; part_num++) {
          if (p[part_num] == 0) continue;
          if (is_proton->at(part_num) == is_pip->at(part_num)) continue;

          theta = physics::theta_calc(cz[part_num]);
          phi = physics::phi_calc(cx[part_num], cy[part_num]);
          sector = physics::get_sector(phi);
          hists->Fill_hadron_fid(theta, phi, sector, id[part_num]);

          hists->Fill_Mass(m[part_num]);
          TLorentzVector Particle = physics::fourVec(p[part_num], cx[part_num], cy[part_num], cz[part_num],
    id[part_num]);

          hists->MomVsBeta_Fill(Particle.E(), p[part_num], b[part_num]);
          if (q[part_num] == 1) {
            hists->MomVsBeta_Fill_pos(p[part_num], b[part_num]);
            if (is_proton->at(part_num) && (id[part_num] == PROTON)) {
              num_of_proton++;
              hists->Fill_proton_WQ2(W, Q2);
              hists->Fill_proton_ID_P(p[part_num], b[part_num]);
            } else if (is_pip->at(part_num) && (id[part_num] == PIP)) {
              num_of_pis++;
              hists->Fill_pion_WQ2(W, Q2);
              hists->Fill_Pi_ID_P(p[part_num], b[part_num]);
              TLorentzVector gamma_mu = (*e_mu - e_mu_prime);
              if (first_run) {
                MM_neutron->Set_PxPyPz(p[part_num] * cx[part_num], p[part_num] * cy[part_num],
                                       p[part_num] * cz[part_num]);
                MM_neutron->missing_mass(gamma_mu);
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

        if (num_of_pis == 1 && gpart < 3) hists->Fill_single_pi_WQ2(W, Q2);

        if (num_of_pis == 1 && num_of_proton == 0 && gpart < 3) {
          hists->Fill_Missing_Mass(MM_neutron->Get_MM());
          hists->Fill_Missing_Mass_square(MM_neutron->Get_MM2());
        }
        if (num_of_proton == 1) hists->Fill_single_proton_WQ2(W, Q2);
      }
      // std::cout << "End of loop " << current_event << std::endl;
    }
    */
  chain->Reset();  // delete Tree object
  delete chain;
}
