/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/
#include "datahandeler.hpp"

DataHandeler::DataHandeler() = default;
DataHandeler::DataHandeler(const std::vector<std::string>& fin, const std::shared_ptr<Histogram>& hists)
    : _input_files(fin), _hists(hists) {}
DataHandeler::~DataHandeler() = default;

void DataHandeler::loadbar(long x, long n) {
  int w = 50;
  if (x >= n) return;

  double ratio = x / (double)n;
  int c = ratio * w;

  std::cout << BLUE << " [";
  for (int x = 0; x < c; x++) std::cout << GREEN << "=" << DEF;
  std::cout << GREEN << ">" << DEF;
  for (int x = c; x < w; x++) std::cout << " ";
  std::cout << BLUE << (int)(ratio * 100) << "%]\r" << DEF << std::flush;
}

void DataHandeler::setLoadBar(bool load) { _loadbar = load; }

int DataHandeler::Run() {
  if (_hists == nullptr) return 0;
  size_t total = 0;
  int current_event = 0;
  int _ci = 0;
  int num_of_events = 0;
#pragma omp parallel private(current_event) num_threads(NUM_THREADS)
  {
#pragma omp critical
    {
      for (_ci = 0; _ci < NUM_THREADS; _ci++) {
        _chain[_ci] = new TChain("h10");
        for (auto& f : _input_files) {
          _chain[_ci]->Add(f.c_str());
        }
        _data[_ci] = std::make_shared<Branches>(_chain[_ci]);
      }

      num_of_events = (int)_chain[omp_get_thread_num()]->GetEntries();
    }
#pragma omp for reduction(+ : total)
    for (current_event = 0; current_event < num_of_events; current_event++) {
      total += DataHandeler::Run(current_event, omp_get_thread_num());
    }
  }
  return total;
}

int DataHandeler::Run(int current_event, int thread) {
#pragma omp critical
  _chain[thread]->GetEntry(current_event);
  auto check = std::make_unique<Cuts>(_data[thread]);
  // if (_data[thread]->ec_eo(0) < 0.01) return 0;
  if (!check->isElecctron()) return 0;

  _hists->EC_inout(_data[thread]->ec_ei(0), _data[thread]->ec_eo(0));
  _hists->EC_fill(_data[thread]->etot(0), _data[thread]->p(0));
  _hists->TM_Fill(_data[thread]->p(0), physics::theta_calc(_data[thread]->cz(0)));
  double theta_cc =
      TMath::ACos(TMath::Abs(_data[thread]->p(0) * _data[thread]->cz(0)) / TMath::Abs(_data[thread]->p(0))) / D2R;
  _hists->CC_fill(_data[thread]->cc_sect(0), (_data[thread]->cc_segm(0) % 1000) / 10,
                  _data[thread]->cc_segm(0) / 1000 - 1, _data[thread]->nphe(0), theta_cc);

  auto event = std::make_unique<Reaction>(_data[thread]);

  _hists->Fill_E_Prime_fid(event->e_mu_prime());
  _hists->Fill_E_Prime(event->e_mu_prime());

  if (getenv("CUTS") != nullptr && atoi(getenv("CUTS")) == true) {
    CUTS = (check->isElecctron() && check->Fid_cut() && check->Beam_cut());
  } else {
    CUTS = check->isElecctron();
  }

  // if (CUTS) {
  if (check->isElecctron()) {
    _hists->Fill_Beam_Position(_data[thread]->dc_vx(0), _data[thread]->dc_vy(0), _data[thread]->dc_vz(0));

    auto dt = std::make_unique<Delta_T>(_data[thread]->sc_t(0), _data[thread]->sc_r(0));

    dt->delta_t_hists(_hists, _data[thread]);
    std::vector<double> dt_proton = dt->delta_t_array(MASS_P, _data[thread]);
    std::vector<double> dt_pi = dt->delta_t_array(MASS_PIP, _data[thread]);

    float theta = physics::theta_calc(_data[thread]->cz(0));
    float phi = physics::phi_calc(_data[thread]->cx(0), _data[thread]->cy(0));
    int sector = _data[thread]->dc_sect(0);
    _hists->Fill_electron_fid(theta, phi, sector);

    // auto photon_flux = std::make_unique<PhotonFlux>(event->e_mu(), event->e_mu_prime());
    //_hists->Photon_flux_Fill(photon_flux->GetVirtualPhotonFlux());

    _hists->WvsQ2_Fill(event->W(), event->Q2(), _data[thread]->ec_sect(0));

    for (int part_num = 1; part_num < _data[thread]->gpart(); part_num++) {
      theta = physics::theta_calc(_data[thread]->cz(part_num));
      phi = physics::phi_calc(_data[thread]->cx(part_num), _data[thread]->cy(part_num));
      sector = _data[thread]->dc_sect(part_num);

      _hists->delta_t_sec_pad(_data[thread]->p(part_num), _data[thread]->q(part_num), dt->Get_dt_P(), dt->Get_dt_Pi(),
                              dt->Get_dt_E(), _data[thread]->sc_sect(part_num), _data[thread]->sc_pd(part_num));

      _hists->Fill_Target_Vertex(_data[thread]->vx(part_num), _data[thread]->vy(part_num), _data[thread]->vz(part_num));
      _hists->MomVsBeta_Fill(_data[thread]->p(part_num), _data[thread]->b(part_num));

      if (_data[thread]->q(part_num) == POSITIVE) {
        _hists->MomVsBeta_Fill_pos(_data[thread]->p(part_num), _data[thread]->b(part_num));
        if (check->dt_Pip_cut(dt_pi.at(part_num), _data[thread]->p(part_num)) &&
            check->dt_P_cut(dt_proton.at(part_num), _data[thread]->p(part_num)))
          _hists->Fill_proton_Pi_ID_P(_data[thread]->p(part_num), _data[thread]->b(part_num));

        if (check->dt_Pip_cut(dt_pi.at(part_num), _data[thread]->p(part_num))) {
          event->SetPip(part_num);
          _hists->Fill_hadron_fid(theta, phi, sector, PIP);
          _hists->Fill_pion_WQ2(event->W(), event->Q2());
          _hists->Fill_Pi_ID_P(_data[thread]->p(part_num), _data[thread]->b(part_num));
        } else if (check->dt_P_cut(dt_proton.at(part_num), _data[thread]->p(part_num))) {
          event->SetProton(part_num);
          _hists->Fill_hadron_fid(theta, phi, sector, PROTON);
          _hists->Fill_proton_WQ2(event->W(), event->Q2());
          _hists->Fill_proton_ID_P(_data[thread]->p(part_num), _data[thread]->b(part_num));
        } else
          event->SetOther(part_num);

      } else if (_data[thread]->q(part_num) == NEGATIVE) {
        _hists->MomVsBeta_Fill_neg(_data[thread]->p(part_num), _data[thread]->b(part_num));
        if (check->dt_Pip_cut(dt_pi.at(part_num), _data[thread]->p(part_num))) {
          _hists->Fill_hadron_fid(theta, phi, sector, PIM);
          event->SetPim(part_num);
        } else
          event->SetOther(part_num);
      } else if (_data[thread]->q(part_num) == 0) {
        _hists->MomVsBeta_Fill_neutral(_data[thread]->p(part_num), _data[thread]->b(part_num));
        event->SetOther(part_num);
      }
    }

    if (event->SingleP() && event->MM() >= 0.1 && event->MM() <= 0.2) _hists->Fill_P_PI0(event->W(), event->Q2());

    if (event->SinglePip()) {
      _hists->Fill_Missing_Mass(event->MM(), event->MM2());
      _hists->Fill_W_Missing_Mass(event->W(), event->MM(), event->MM2());
    }
    bool mm_cut = true;

    mm_cut &= (event->MM() < 0.987669);
    mm_cut &= (event->MM() > 0.923374);
    if (mm_cut) _hists->Fill_MM_WQ2(event->W(), event->Q2());
    if ((event->SinglePip() || event->NeutronPip()) && mm_cut) {
      _hists->Fill_channel_WQ2(event->W(), event->Q2(), _data[thread]->ec_sect(0), event->e_mu_prime(), event->MM(),
                               event->MM2());
      _hists->Fill_Missing_Mass_strict(event->MM(), event->MM2());
      _hists->EC_cut_fill(_data[thread]->etot(0), _data[thread]->p(0));
      _hists->Fill_E_Prime_channel(event->e_mu_prime());
    }
    if ((event->SinglePip() || event->NeutronPip()))
      _hists->Fill_NeutronPip_WQ2(event->W(), event->Q2(), event->MM(), event->MM2());
    if (event->SingleP()) {
      _hists->Fill_single_proton_WQ2(event->W(), event->Q2());
      _hists->Fill_Missing_Mass_pi0(event->MM(), event->MM2());
    }
    if (event->TwoPion()) _hists->Fill_Missing_Mass_twoPi(event->MM(), event->MM2());
  }

  return 1;
}

mcHandeler::mcHandeler(const std::vector<std::string>& fin, const std::shared_ptr<mcHistogram>& hists) {
  _input_files = fin;
  _hists = hists;
  _mc_hists = hists;
}

int mcHandeler::Run() {
  size_t total = 0;
  int _ci = 0;
  for (auto& _c : _chain) {
    _c = new TChain("h10");
    for (auto& f : _input_files) {
      _c->Add(f.c_str());
    }
    _data[_ci++] = std::make_shared<Branches>(_c, true);
  }

  int num_of_events = (int)_chain[0]->GetEntries();
  int current_event = 0;
#pragma omp parallel for private(current_event)
  for (current_event = 0; current_event < num_of_events; current_event++) {
#pragma omp critical
    if (current_event % 500000 == 0 && _loadbar) DataHandeler::loadbar(current_event + 1, num_of_events);
    total += DataHandeler::Run(current_event, omp_get_thread_num());
    total += mcHandeler::Run(current_event, omp_get_thread_num());
  }
  //_chain->Reset();
  return total;
}

int mcHandeler::Run(int current_event, int thread) {
#pragma omp critical
  _chain[thread]->GetEntry(current_event);
  auto check = std::make_unique<Cuts>(_data[thread]);
  // if (!check->isElecctron()) return 0;

  auto mc_event = std::make_unique<Reaction>(_data[thread], true);
  _mc_hists->Fill_P(_data[thread]);
  _mc_hists->Fill_WQ2_MC(mc_event->W(), mc_event->Q2());

  for (int part_num = 1; part_num < _data[thread]->gpart(); part_num++) {
    if (_data[thread]->pidpart(part_num) == PIP) mc_event->SetPip(part_num);
  }
  if (mc_event->SinglePip()) _mc_hists->Fill_Missing_Mass(mc_event->MM(), mc_event->MM2());

  return 1;
}
