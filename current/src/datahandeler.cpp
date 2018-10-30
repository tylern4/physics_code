/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#include "datahandeler.hpp"

DataHandeler::DataHandeler() {
  if (getenv("BEAM_E") != NULL) {
    BEAM_ENERGY = atof(getenv("BEAM_E"));
    std::cout << RED << "Beam energy set to: " << BEAM_ENERGY << DEF << std::endl;
  } else {
    BEAM_ENERGY = E1D_E0;
  }
}
DataHandeler::~DataHandeler() {}

void DataHandeler::Run(std::vector<std::string> fin, Histogram *hists) {
  for (auto f : fin) Run(f, hists);
}

void DataHandeler::Run(std::string fin, Histogram *hists) {
  TLorentzVector *e_mu = new TLorentzVector(0.0, 0.0, sqrt(BEAM_ENERGY * BEAM_ENERGY - MASS_E * MASS_E), BEAM_ENERGY);

  MissingMass *MM_neutron = new MissingMass(MASS_P, 0.0);
  MissingMass *MM_pi0 = new MissingMass(MASS_P, 0.0);
  MissingMass *MM_from2pi = new MissingMass(MASS_P, 0.0);
  // End declrare variables

  // Load chain from branch h10
  bool cuts, electron_cuts;
  int num_of_events;
  int total_events;
  int num_of_pips, num_of_pims;
  int num_of_proton, neg;
  double W, Q2;
  double e_E;
  double theta;
  double phi;
  int sector;
  // bool first_run = true;
  TChain *chain = new TChain("h10");
  chain->Add(fin.c_str());
  Branches *data = new Branches(chain);
  num_of_events = (int)chain->GetEntries();
  int current_event = 0;
  for (current_event = 0; current_event < num_of_events; current_event++) {
    chain->GetEntry(current_event);
    Cuts *check = new Cuts();
    theta = physics::theta_calc(data->cz(0));
    phi = physics::phi_calc(data->cx(0), data->cy(0));
    sector = physics::get_sector(phi);
    check->Set_elec_fid(theta, phi, sector);
    check->Set_BeamPosition(data->dc_vx(data->dc(0) - 1), data->dc_vy(data->dc(0) - 1), data->dc_vz(data->dc(0) - 1));

    // Setup scattered electron 4 vector
    TLorentzVector e_mu_prime = physics::fourVec(data->p(0), data->cx(0), data->cy(0), data->cz(0), MASS_E);
    W = physics::W_calc(*e_mu, e_mu_prime);
    Q2 = physics::Q2_calc(*e_mu, e_mu_prime);

    hists->Fill_E_Prime(e_mu_prime);

    if (check->Fid_cut()) {
      hists->EC_fill(data->etot(data->ec(0) - 1), data->p(0));
      hists->Fill_E_Prime_fid(e_mu_prime);
      hists->TM_Fill(data->p(0), physics::theta_calc(data->cz(0)));
    }

    if (getenv("CUTS") != NULL && atoi(getenv("CUTS")) == true) {
      CUTS = (check->Fid_cut() && check->Beam_cut());
    } else {
      CUTS = true;
    }

    if (CUTS) {
      int cc_sector = data->cc_sect(data->cc(0) - 1);
      int cc_segment = (data->cc_segm(0) % 1000) / 10;
      int cc_pmt = data->cc_segm(0) / 1000 - 1;
      int cc_nphe = data->nphe(data->cc(0) - 1);
      double theta_cc = TMath::ACos(TMath::Abs(data->p(0) * data->cz(0)) / TMath::Abs(data->p(0)));

      theta_cc = theta_cc / D2R;
      hists->CC_fill(cc_sector, cc_segment, cc_pmt, cc_nphe, theta_cc);
      hists->Fill_Beam_Position(data->dc_vx(data->dc(0) - 1), data->dc_vy(data->dc(0) - 1),
                                data->dc_vz(data->dc(0) - 1));

      // Set the vertex time (time of electron hit)
      Delta_T *dt = new Delta_T(data->sc_t(data->sc(0) - 1), data->sc_r(data->sc(0) - 1));
      dt->delta_t_hists(hists, data);
      std::vector<double> dt_proton = dt->delta_t_array(MASS_P, data);
      std::vector<double> dt_pi = dt->delta_t_array(MASS_PIP, data);
      delete dt;
      hists->Fill_electron_fid(theta, phi, sector);

      e_E = e_mu_prime.E();
      PhotonFlux *photon_flux = new PhotonFlux(*e_mu, e_mu_prime);
      hists->Photon_flux_Fill(photon_flux->GetVirtualPhotonFlux());
      delete photon_flux;
      TLorentzVector gamma_mu = (*e_mu - e_mu_prime);
      hists->WvsQ2_Fill(W, Q2);
      num_of_proton = num_of_pips = num_of_pims = neg = 0;
      TLorentzVector particle;
      for (int part_num = 1; part_num < data->gpart(); part_num++) {
        theta = physics::theta_calc(data->cz(part_num));
        phi = physics::phi_calc(data->cx(part_num), data->cy(part_num));
        sector = physics::get_sector(phi);

        if (data->q(part_num) == POSITIVE && check->dt_Pip_cut(dt_pi.at(part_num), data->p(part_num))) {
          hists->Fill_hadron_fid(theta, phi, sector, PIP);
          particle =
              physics::fourVec(data->p(part_num), data->cx(part_num), data->cy(part_num), data->cz(part_num), MASS_PIP);
        } else if (data->q(part_num) == POSITIVE && check->dt_P_cut(dt_proton.at(part_num), data->p(part_num))) {
          hists->Fill_hadron_fid(theta, phi, sector, PROTON);
          particle =
              physics::fourVec(data->p(part_num), data->cx(part_num), data->cy(part_num), data->cz(part_num), MASS_P);
        } else if (data->q(part_num) == NEGATIVE && check->dt_Pip_cut(dt_pi.at(part_num), data->p(part_num))) {
          hists->Fill_hadron_fid(theta, phi, sector, PIM);
          particle =
              physics::fourVec(data->p(part_num), data->cx(part_num), data->cy(part_num), data->cz(part_num), MASS_PIM);
        } else {
          continue;
        }

        hists->Fill_Target_Vertex(data->vx(part_num), data->vy(part_num), data->vz(part_num));
        hists->MomVsBeta_Fill(particle.E(), data->p(part_num), data->b(part_num));
        if (data->q(part_num) == POSITIVE) {
          hists->MomVsBeta_Fill_pos(data->p(part_num), data->b(part_num));
          if (check->dt_Pip_cut(dt_pi.at(part_num), data->p(part_num)) &&
              check->dt_P_cut(dt_proton.at(part_num), data->p(part_num)))
            hists->Fill_proton_Pi_ID_P(data->p(part_num), data->b(part_num));

          if (check->dt_Pip_cut(dt_pi.at(part_num), data->p(part_num))) {
            num_of_pips++;
            hists->Fill_pion_WQ2(W, Q2);
            hists->Fill_Pi_ID_P(data->p(part_num), data->b(part_num));
            MM_neutron->Set_4Vec(particle);
            MM_neutron->missing_mass(gamma_mu);
          } else if (check->dt_P_cut(dt_proton.at(part_num), data->p(part_num))) {
            particle =
                physics::fourVec(data->p(part_num), data->cx(part_num), data->cy(part_num), data->cz(part_num), PROTON);
            num_of_proton++;
            hists->Fill_proton_WQ2(W, Q2);
            hists->Fill_proton_ID_P(data->p(part_num), data->b(part_num));
            MM_pi0->Set_4Vec(particle);
            MM_pi0->missing_mass(gamma_mu);
          }
        } else if (data->q(part_num) == NEGATIVE) {
          neg++;
          hists->MomVsBeta_Fill_neg(data->p(part_num), data->b(part_num));
          if (check->dt_Pip_cut(dt_pi.at(part_num), data->p(part_num))) {
            particle =
                physics::fourVec(data->p(part_num), data->cx(part_num), data->cy(part_num), data->cz(part_num), PIM);
            num_of_pims++;
            MM_from2pi->Set_4Vec(particle);
            MM_from2pi->missing_mass(gamma_mu);
          }
        }
      }

      bool mm_cut = true;
      mm_cut &= (MM_neutron->Get_MM() < 1.1);
      mm_cut &= (MM_neutron->Get_MM() > 0.8);

      if (num_of_proton == 1 && num_of_pips == 0 && num_of_pims == 0 && MM_pi0->Get_MM() >= 0.1 &&
          MM_pi0->Get_MM() <= 0.2)
        hists->Fill_P_PI0(W, Q2);

      if (mm_cut) hists->Fill_MM_WQ2(W, Q2);
      if (num_of_pips == 1 && num_of_proton == 0 && num_of_pims == 0) {
        hists->Fill_Missing_Mass(MM_neutron);
      }
      if (num_of_pips == 1 && num_of_proton == 0 && num_of_pims == 0 && mm_cut && neg == 0) {
        hists->Fill_channel_WQ2(W, Q2, e_mu_prime, MM_neutron->Get_MM(), MM_neutron->Get_MM2(), sector);
        hists->Fill_Missing_Mass_strict(MM_neutron);
        hists->EC_cut_fill(data->etot(data->ec(0) - 1), data->p(0));
        hists->Fill_E_Prime_channel(e_mu_prime);
      }
      if (num_of_pips == 1 && num_of_proton == 0 && num_of_pims == 0) hists->Fill_single_pi_WQ2(W, Q2);
      if (num_of_pips == 0 && num_of_proton == 1 && num_of_pims == 0) hists->Fill_single_proton_WQ2(W, Q2);
      if (num_of_proton == 1) hists->Fill_Missing_Mass_pi0(MM_pi0);
      if (num_of_pips == 1 && num_of_pims == 1) hists->Fill_Missing_Mass_twoPi(MM_from2pi);
    }
    delete check;
  }
  chain->Reset();  // delete Tree object
}

void DataHandeler::loadbar(long x, long n) {
  int w = 50;
  if (x <= n) return;
  if ((x != n) && (x % (n / 100 + 1) != 0)) return;

  double ratio = x / (double)n;
  int c = ratio * w;

  std::cout << BLUE << " [";
  for (int x = 0; x < c; x++) std::cout << GREEN << "=" << DEF;
  std::cout << GREEN << ">" << DEF;
  for (int x = c; x < w; x++) std::cout << " ";
  std::cout << BLUE << (int)(ratio * 100) << "%]\r" << DEF << std::flush;
}
