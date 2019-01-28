/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#include "datahandeler.hpp"

DataHandeler::DataHandeler() {}
DataHandeler::~DataHandeler() {}

void DataHandeler::Run(std::vector<std::string> fin, Histogram *hists) {
  for (auto f : fin) Run(f, hists);
}

void DataHandeler::Run(std::string fin, Histogram *hists) {
  TChain *chain = new TChain("h10");
  chain->Add(fin.c_str());
  Branches *data = new Branches(chain);
  int num_of_events = (int)chain->GetEntries();

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain->GetEntry(current_event);
    auto check = std::make_unique<Cuts>(data);
    if (!check->isElecctron()) continue;
    if (data->ec_ei(0) < 0.01 || data->ec_eo(0) < 0.01) continue;

    hists->EC_inout(data->ec_ei(0), data->ec_eo(0));
    hists->EC_fill(data->etot(0), data->p(0));
    hists->TM_Fill(data->p(0), physics::theta_calc(data->cz(0)));
    double theta_cc = TMath::ACos(TMath::Abs(data->p(0) * data->cz(0)) / TMath::Abs(data->p(0))) / D2R;
    hists->CC_fill(data->cc_sect(0), (data->cc_segm(0) % 1000) / 10, data->cc_segm(0) / 1000 - 1, data->nphe(0),
                   theta_cc);

    auto event = std::make_unique<Reaction>(data);

    hists->Fill_E_Prime_fid(event->e_mu_prime());
    hists->Fill_E_Prime(event->e_mu_prime());

    if (getenv("CUTS") != NULL && atoi(getenv("CUTS")) == true) {
      CUTS = (check->isElecctron() && check->Fid_cut() && check->Beam_cut());
    } else {
      CUTS = check->isElecctron();
    }

    if (check->isElecctron()) {
      hists->Fill_Beam_Position(data->dc_vx(0), data->dc_vy(0), data->dc_vz(0));

      auto dt = std::make_unique<Delta_T>(data->sc_t(0), data->sc_r(0));

      dt->delta_t_hists(hists, data);
      std::vector<double> dt_proton = dt->delta_t_array(MASS_P, data);
      std::vector<double> dt_pi = dt->delta_t_array(MASS_PIP, data);

      float theta = physics::theta_calc(data->cz(0));
      float phi = physics::phi_calc(data->cx(0), data->cy(0));
      int sector = data->dc_sect(0);
      hists->Fill_electron_fid(theta, phi, sector);

      auto photon_flux = std::make_unique<PhotonFlux>(event->e_mu(), event->e_mu_prime());
      hists->Photon_flux_Fill(photon_flux->GetVirtualPhotonFlux());

      hists->WvsQ2_Fill(event->W(), event->Q2(), data->ec_sect(0));

      for (int part_num = 1; part_num < data->gpart(); part_num++) {
        theta = physics::theta_calc(data->cz(part_num));
        phi = physics::phi_calc(data->cx(part_num), data->cy(part_num));
        sector = data->dc_sect(part_num);

        hists->Fill_Target_Vertex(data->vx(part_num), data->vy(part_num), data->vz(part_num));
        hists->MomVsBeta_Fill(data->p(part_num), data->b(part_num));

        if (data->q(part_num) == POSITIVE) {
          hists->MomVsBeta_Fill_pos(data->p(part_num), data->b(part_num));
          if (check->dt_Pip_cut(dt_pi.at(part_num), data->p(part_num)) &&
              check->dt_P_cut(dt_proton.at(part_num), data->p(part_num)))
            hists->Fill_proton_Pi_ID_P(data->p(part_num), data->b(part_num));

          if (check->dt_Pip_cut(dt_pi.at(part_num), data->p(part_num))) {
            event->SetPip(part_num);
            hists->Fill_hadron_fid(theta, phi, sector, PIP);
            hists->Fill_pion_WQ2(event->W(), event->Q2());
            hists->Fill_Pi_ID_P(data->p(part_num), data->b(part_num));
          } else if (check->dt_P_cut(dt_proton.at(part_num), data->p(part_num))) {
            event->SetProton(part_num);
            hists->Fill_hadron_fid(theta, phi, sector, PROTON);
            hists->Fill_proton_WQ2(event->W(), event->Q2());
            hists->Fill_proton_ID_P(data->p(part_num), data->b(part_num));
          } else
            event->SetOther(part_num);

        } else if (data->q(part_num) == NEGATIVE) {
          hists->MomVsBeta_Fill_neg(data->p(part_num), data->b(part_num));
          if (check->dt_Pip_cut(dt_pi.at(part_num), data->p(part_num))) {
            hists->Fill_hadron_fid(theta, phi, sector, PIM);
            event->SetPim(part_num);
          } else
            event->SetOther(part_num);
        } else if (data->q(part_num) == 0) {
          hists->MomVsBeta_Fill_neutral(data->p(part_num), data->b(part_num));
          event->SetOther(part_num);
        }
      }

      if (event->SingleP() && event->MM() >= 0.1 && event->MM() <= 0.2) hists->Fill_P_PI0(event->W(), event->Q2());

      if (event->SinglePip()) {
        hists->Fill_Missing_Mass(event->MM(), event->MM2());
        hists->Fill_W_Missing_Mass(event->W(), event->MM(), event->MM2());
      }
      bool mm_cut = true;
      // mm_cut &= (event->MM() < 1.2);
      // mm_cut &= (event->MM() > 0.8);

      mm_cut &= (event->MM() < 0.987669);
      mm_cut &= (event->MM() > 0.923374);
      if (mm_cut) hists->Fill_MM_WQ2(event->W(), event->Q2());
      if ((event->SinglePip() || event->NeutronPip()) && mm_cut) {
        hists->Fill_channel_WQ2(event->W(), event->Q2(), data->ec_sect(0), event->e_mu_prime(), event->MM(),
                                event->MM2());
        hists->Fill_Missing_Mass_strict(event->MM(), event->MM2());
        hists->EC_cut_fill(data->etot(0), data->p(0));
        hists->Fill_E_Prime_channel(event->e_mu_prime());
      }
      if ((event->SinglePip() || event->NeutronPip())) hists->Fill_NeutronPip_WQ2(event->W(), event->Q2());
      if (event->SingleP()) {
        hists->Fill_single_proton_WQ2(event->W(), event->Q2());
        hists->Fill_Missing_Mass_pi0(event->MM(), event->MM2());
      }
      if (event->TwoPion()) hists->Fill_Missing_Mass_twoPi(event->MM(), event->MM2());
    }
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
