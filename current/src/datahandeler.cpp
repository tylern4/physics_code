/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#include "datahandeler.hpp"

DataHandeler::DataHandeler(std::vector<std::string> fin, std::string RootFile_output) {
  input_files = fin;
  c1 = new TCanvas("c1", "c1", 100, 100);

  double BEAM_ENERGY;
  if (getenv("BEAM_E") != NULL) {
    BEAM_ENERGY = atof(getenv("BEAM_E"));
  } else {
    BEAM_ENERGY = E1D_E0;
  }
  e_mu = new TLorentzVector(0.0, 0.0, sqrt(BEAM_ENERGY * BEAM_ENERGY - MASS_E * MASS_E), BEAM_ENERGY);
  // End declrare variables

  // Open outputfile
  RootOutputFile = new TFile(RootFile_output.c_str(), "RECREATE");
  hists = new Histogram();
  MM_neutron = new MissingMass(MASS_P, 0.0);
  MM_pi0 = new MissingMass(MASS_P, 0.0);
  MM_from2pi = new MissingMass(MASS_P, 0.0);
}

DataHandeler::~DataHandeler() {
  delete MM_neutron;
  std::cout << GREEN << "\nFitting" << DEF << std::endl;
  // Start of cuts
  Fits *MM_neutron_cut = new Fits();
  MM_neutron_cut->Set_min(0.8);
  MM_neutron_cut->Set_max(1.2);
  MM_neutron_cut->FitBreitWigner(hists->Missing_Mass_strict);
  // MM_neutron_cut->FitGaus(hists->Missing_Mass_strict);
  // MM_neutron_cut->Fit2Gaus(hists->Missing_Mass_strict);
  // MM_neutron_cut->FitLandau(hists->Missing_Mass_strict);

  // Header *MM_header = new Header("../src/missing_mass_gaussians.hpp", "MM");
  // MM_header->WriteGaussian("mm", 1, MM_neutron_cut->Get_mean(), MM_neutron_cut->Get_sigma());

  Fits *MissingMassSquare_cut = new Fits();
  MissingMassSquare_cut->Set_min(0.5);
  MissingMassSquare_cut->Set_max(1.1);
  MissingMassSquare_cut->FitBreitWigner(hists->Missing_Mass_square_strict);
  // MissingMassSquare_cut->FitGaus(hists->Missing_Mass_square_strict);
  // MissingMassSquare_cut->Fit2Gaus(hists->Missing_Mass_square_strict);
  // MissingMassSquare_cut->FitLandau(hists->Missing_Mass_square_strict);
  // MM_header->WriteGaussian("mm_square", 1, MissingMassSquare_cut->Get_mean(), MissingMassSquare_cut->Get_sigma());
  // delete MM_header;
  delete MM_neutron_cut;

  //
  // end stuff

  RootOutputFile->cd();
  std::cerr << BOLDBLUE << "EC_Write()" << DEF << std::endl;
  TDirectory *EC_folder = RootOutputFile->mkdir("EC_hists");
  EC_folder->cd();
  hists->EC_Write();
  std::cerr << BOLDBLUE << "EC_slices()" << DEF << std::endl;
  TDirectory *EC_slices = RootOutputFile->mkdir("EC_slices");
  EC_slices->cd();
  hists->EC_slices_Write();
  std::cerr << BOLDBLUE << "Beam_Position()" << DEF << std::endl;
  TDirectory *Beam_Folder = RootOutputFile->mkdir("Beam Position");
  Beam_Folder->cd();
  hists->Beam_Position_Write();
  hists->Target_Vertex_Write();
  std::cerr << BOLDBLUE << "WvsQ2()" << DEF << std::endl;
  TDirectory *WvsQ2_folder = RootOutputFile->mkdir("W vs Q2");
  WvsQ2_folder->cd();
  hists->WvsQ2_Write();
  TDirectory *W_Q2_binned = RootOutputFile->mkdir("W_Q2_binned");
  W_Q2_binned->cd();
  hists->WvsQ2_binned_Write();
  std::cerr << BOLDBLUE << "MomVsBeta_Fill()" << DEF << std::endl;
  TDirectory *MomVsBeta_folder = RootOutputFile->mkdir("Momentum vs beta");
  MomVsBeta_folder->cd();
  hists->MomVsBeta_Write();
  std::cerr << BOLDBLUE << "Write_Missing_Mass()" << DEF << std::endl;
  // Missing Mass Write
  TDirectory *MissMass = RootOutputFile->mkdir("Missing_Mass");
  MissMass->cd();
  hists->Write_Missing_Mass();
  std::cerr << BOLDBLUE << "delta_t_Write()" << DEF << std::endl;
  // Delta T Write
  TDirectory *DeltaT = RootOutputFile->mkdir("Delta_T");
  DeltaT->cd();
  hists->delta_t_Write();
  std::cerr << BOLDBLUE << "delta_t_slices_Write()" << DEF << std::endl;
  TDirectory *DeltaT_slices = RootOutputFile->mkdir("Delta_T_slices");
  DeltaT_slices->cd();
  hists->delta_t_slices_Write();
  std::cerr << BOLDBLUE << "delta_t_sec_pad_Write()" << DEF << std::endl;
  TDirectory *DeltaT_sec_pad = RootOutputFile->mkdir("Delta_T_sec_pad");
  DeltaT_sec_pad->cd();
  hists->delta_t_sec_pad_Write();
  std::cerr << BOLDBLUE << "delta_T_canvas()" << DEF << std::endl;
  TDirectory *Delta_T_canvases = RootOutputFile->mkdir("Delta_T_canvases");
  Delta_T_canvases->cd();
  hists->delta_T_canvas();
  std::cerr << BOLDBLUE << "Theta_CC_Write()" << DEF << std::endl;
  TDirectory *Theta_CC_hists = RootOutputFile->mkdir("Theta_CC_hists");
  Theta_CC_hists->cd();
  hists->Theta_CC_Write();
  std::cerr << BOLDBLUE << "CC_Write()" << DEF << std::endl;
  TDirectory *CC_hists = RootOutputFile->mkdir("CC_hists");
  CC_hists->cd();
  hists->CC_Write();
  std::cerr << BOLDBLUE << "CC_canvas()" << DEF << std::endl;
  TDirectory *CC_canvases = RootOutputFile->mkdir("CC_canvases");
  CC_canvases->cd();
  hists->CC_canvas();
  std::cerr << BOLDBLUE << "Fid_Write()" << DEF << std::endl;
  TDirectory *Fid_cuts = RootOutputFile->mkdir("Fid_cuts");
  Fid_cuts->cd();
  hists->Fid_Write();
  std::cerr << BOLDBLUE << "fid_canvas()" << DEF << std::endl;
  TDirectory *Fid_canvas = RootOutputFile->mkdir("Fid_canvas");
  Fid_canvas->cd();
  hists->fid_canvas();
  std::cerr << BOLDBLUE << "Done!!!" << DEF << std::endl;
  delete hists;
  // RootOutputFile->Write();
  RootOutputFile->Close();
  // cut_outputs.close();
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
  std::thread *fh_thread[size];

  for (i = 0; i < size; i++) {
    loadbar(i, size - 1);
#ifndef __THREAD__
    file_handeler(input_files.at(i));
#else
    try {
      fh_thread[i] = new std::thread(std::mem_fn(&DataHandeler::file_handeler), this, input_files.at(i));
      fh_thread[i]->join();
    } catch (const std::exception &e) {
      std::cerr << RED << "Error:\t" << e.what() << std::endl;
      std::cerr << CYAN << "Bad File: \t" << input_files.at(i) << DEF << std::endl;
    }
#endif
  }
}

void DataHandeler::file_handeler(std::string fin) {
  // Load chain from branch h10
  bool cuts, electron_cuts;
  int num_of_events;
  int total_events;
  int num_of_pips, num_of_pims;
  int num_of_proton;
  double e_E;
  double theta;
  double phi;
  int sector;
  // bool first_run = true;

  TChain *chain = new TChain("h10");
  chain->Add(fin.c_str());

  getBranches(chain);
  // if (!first_run) getMorebranchs(chain);
  num_of_events = (int)chain->GetEntries();

  int current_event = 0;
  for (current_event = 0; current_event < num_of_events; current_event++) {
    chain->GetEntry(current_event);
    // if (gpart >= 3) continue;
    // if (p[0] < 1.0) continue;
    // if (abs((double)dc_vz[dc[0] - 1]) > 2) continue;
    Cuts *check = new Cuts();

    // electron cuts
    check->Set_charge((int)q[0]);
    check->Set_ec_cut(ec[0] > 0);        // ``` ``` ``` ec
    check->Set_electron_id((int)id[0]);  // First particle is electron
    check->Set_gpart((int)gpart);        // Number of good particles is greater than 0
    check->Set_cc_cut((int)cc[0] > 0);
    check->Set_stat_cut((int)stat[0] > 0);  // First Particle hit stat
    check->Set_sc_cut((int)sc[0] > 0);
    check->Set_dc_cut((int)dc[0] > 0);
    check->Set_dc_stat_cut((int)dc_stat[dc[0] - 1] > 0);
    check->Set_p((double)p[0]);
    if (check->isElecctron()) hists->EC_fill(etot[ec[0] - 1], p[0]);
    if (check->isElecctron()) hists->TM_Fill(p[0], physics::theta_calc(cz[0]));
    check->Set_Sf((double)etot[ec[0] - 1] / (double)p[0]);
    check->Set_num_phe((int)nphe[cc[0] - 1]);
    // Beam position cut
    check->Set_BeamPosition((double)dc_vx[dc[0] - 1], (double)dc_vy[dc[0] - 1], (double)dc_vz[dc[0] - 1]);
    if (check->isStrictElecctron()) {
      int cc_sector = cc_sect[cc[0] - 1];
      int cc_segment = (cc_segm[0] % 1000) / 10;
      int cc_pmt = cc_segm[0] / 1000 - 1;
      int cc_nphe = nphe[cc[0] - 1];
      double theta_cc = TMath::ACos(TMath::Abs(p[0] * cz[0]) / TMath::Abs(p[0]));

      theta_cc = theta_cc / D2R;
      hists->CC_fill(cc_sector, cc_segment, cc_pmt, cc_nphe, theta_cc);

      hists->EC_cut_fill(etot[ec[0] - 1], p[0]);
      hists->Fill_Beam_Position((double)dc_vx[dc[0] - 1], (double)dc_vy[dc[0] - 1], (double)dc_vz[dc[0] - 1]);

      // Setup scattered electron 4 vector
      TLorentzVector e_mu_prime = physics::fourVec(p[0], cx[0], cy[0], cz[0], MASS_E);
      // Set the vertex time (time of electron hit)
      Delta_T *dt = new Delta_T(sc_t[sc[0] - 1], sc_r[sc[0] - 1]);
      dt->delta_t_hists(hists);
      std::vector<double> dt_proton = dt->delta_t_array(MASS_P, gpart);
      std::vector<double> dt_pi = dt->delta_t_array(MASS_PIP, gpart);
      delete dt;
      theta = physics::theta_calc(cz[0]);
      phi = physics::phi_calc(cx[0], cy[0]);
      sector = physics::get_sector(phi);
      hists->Fill_electron_fid(theta, phi, sector);
      W = physics::W_calc(*e_mu, e_mu_prime);
      Q2 = physics::Q2_calc(*e_mu, e_mu_prime);
      e_E = e_mu_prime.E();
      PhotonFlux *photon_flux = new PhotonFlux(*e_mu, e_mu_prime);
      hists->Photon_flux_Fill(photon_flux->GetVirtualPhotonFlux());
      delete photon_flux;
      TLorentzVector gamma_mu = (*e_mu - e_mu_prime);
      hists->WvsQ2_Fill(W, Q2);
      num_of_proton = num_of_pips = num_of_pims = 0;
      for (int part_num = 1; part_num < gpart; part_num++) {
        if (gpart > 3) continue;
        if (p[part_num] == 0) continue;
        PID = -99;
        if (q[part_num] == POSITIVE) {
          if (check->dt_P_cut(dt_proton.at(part_num), p[part_num]))
            PID = PROTON;
          else if (check->dt_Pip_cut(dt_proton.at(part_num), p[part_num]))
            PID = PIP;
        } else if (q[part_num] == NEGATIVE) {
          if (check->dt_Pip_cut(dt_proton.at(part_num), p[part_num])) PID = PIM;
        }

        TLorentzVector Particle = physics::fourVec(p[part_num], cx[part_num], cy[part_num], cz[part_num], PID);
        hists->Fill_Target_Vertex((double)vx[part_num], (double)vy[part_num], (double)vz[part_num]);

        theta = physics::theta_calc(cz[part_num]);
        phi = physics::phi_calc(cx[part_num], cy[part_num]);

        sector = physics::get_sector(phi);
        hists->Fill_hadron_fid(theta, phi, sector, PID);
        hists->MomVsBeta_Fill(Particle.E(), p[part_num], b[part_num]);
        if (q[part_num] == POSITIVE) {
          hists->MomVsBeta_Fill_pos(p[part_num], b[part_num]);
          if (check->dt_P_cut(dt_proton.at(part_num), p[part_num])) {
            num_of_proton++;
            hists->Fill_proton_WQ2(W, Q2);
            hists->Fill_proton_ID_P(p[part_num], b[part_num]);
            MM_pi0->Set_4Vec(Particle);
            MM_pi0->missing_mass(gamma_mu);

          } else if (check->dt_Pip_cut(dt_pi.at(part_num), p[part_num])) {
            num_of_pips++;
            hists->Fill_pion_WQ2(W, Q2);
            hists->Fill_Pi_ID_P(p[part_num], b[part_num]);
            MM_neutron->Set_4Vec(Particle);
            MM_neutron->missing_mass(gamma_mu);
          }

          if (check->dt_Pip_cut(dt_pi.at(part_num), p[part_num]) &&
              check->dt_P_cut(dt_proton.at(part_num), p[part_num])) {
            hists->Fill_proton_Pi_ID_P(p[part_num], b[part_num]);
          }
        } else if (q[part_num] == NEGATIVE) {
          if (check->dt_Pip_cut(dt_pi.at(part_num), p[part_num])) {
            num_of_pims++;
            MM_from2pi->missing_mass(gamma_mu);
            MM_from2pi->Set_4Vec(Particle);
          }
          hists->MomVsBeta_Fill_neg(p[part_num], b[part_num]);
        }
      }

      bool mm_cut = true;
      mm_cut &= (MM_neutron->Get_MM() < 1.1);
      mm_cut &= (MM_neutron->Get_MM() > 0.8);

      if (num_of_pips == 1 && gpart == 2) hists->Fill_single_pi_WQ2(W, Q2);
      if (num_of_proton == 1 && gpart == 2) hists->Fill_single_proton_WQ2(W, Q2);
      if (num_of_pips == 1) hists->Fill_Missing_Mass(MM_neutron);
      if (num_of_pips == 1 && mm_cut && num_of_proton == 0 && gpart < 3) {
        hists->Fill_channel_WQ2(W, Q2, e_mu_prime.E(), physics::xb_calc(Q2, e_mu_prime.E()));
        hists->Fill_Missing_Mass_strict(MM_neutron);
      }
      if (num_of_pips == 2) hists->Fill_Missing_Mass_twoPi(MM_from2pi);
      if (num_of_proton == 1) hists->Fill_Missing_Mass_pi0(MM_pi0);
    }
    delete check;
  }
  chain->Reset();  // delete Tree object
}

void DataHandeler::make_events() {
  int num_of_events;
  TChain *chain = new TChain("h10");
  for (auto fin : input_files) chain->Add(fin.c_str());

  getBranches(chain);
  num_of_events = (int)chain->GetEntries();

  int current_event = 0;
  for (current_event = 0; current_event < num_of_events; current_event++) {
    chain->GetEntry(current_event);
    Cuts *check = new Cuts();

    // electron cuts
    check->Set_charge((int)q[0]);
    check->Set_ec_cut(ec[0] > 0);        // ``` ``` ``` ec
    check->Set_electron_id((int)id[0]);  // First particle is electron
    check->Set_gpart((int)gpart);        // Number of good particles is greater than 0
    check->Set_cc_cut((int)cc[0] > 0);
    check->Set_stat_cut((int)stat[0] > 0);  // First Particle hit stat
    check->Set_sc_cut((int)sc[0] > 0);
    check->Set_dc_cut((int)dc[0] > 0);
    check->Set_dc_stat_cut((int)dc_stat[dc[0] - 1] > 0);
    check->Set_p((double)p[0]);
    check->Set_Sf((double)etot[ec[0] - 1] / (double)p[0]);
    check->Set_num_phe((int)nphe[cc[0] - 1]);
    // Beam position cut
    check->Set_BeamPosition((double)dc_vx[dc[0] - 1], (double)dc_vy[dc[0] - 1], (double)dc_vz[dc[0] - 1]);

    if (check->isStrictElecctron()) {
      Delta_T *dt = new Delta_T(sc_t[sc[0] - 1], sc_r[sc[0] - 1]);
      dt->delta_t_hists(hists);
      std::vector<double> dt_proton = dt->delta_t_array(MASS_P, gpart);
      std::vector<double> dt_pi = dt->delta_t_array(MASS_PIP, gpart);
      delete dt;

      Particle *elec = new Particle(p[0], cx[0], cy[0], cz[0], ELECTRON);
      Event *events = new Event(*elec);

      for (int part_num = 1; part_num < gpart; part_num++) {
        if (q[part_num] == POSITIVE) {
          if (check->dt_P_cut(dt_proton.at(part_num), p[part_num])) {
            Particle part(p[part_num], cx[part_num], cy[part_num], cz[part_num], PROTON);
            events->Add_Part(part);
          } else if (check->dt_Pip_cut(dt_pi.at(part_num), p[part_num])) {
            Particle part(p[part_num], cx[part_num], cy[part_num], cz[part_num], PIP);
            events->Add_Part(part);
          }
        } else if (q[part_num] == NEGATIVE && check->dt_Pip_cut(dt_pi.at(part_num), p[part_num])) {
          Particle part(p[part_num], cx[part_num], cy[part_num], cz[part_num], PIM);
          events->Add_Part(part);
        }
      }
      events->PrintSigniture();
    }
    delete check;
  }
  chain->Reset();  // delete Tree object
}
