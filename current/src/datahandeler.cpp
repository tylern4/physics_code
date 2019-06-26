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

  _chain = new TChain("h10");
  for (auto& f : _input_files) _chain->Add(f.c_str());
  _data = std::make_shared<Branches>(_chain);
  num_of_events = (int)_chain->GetEntries();

  for (current_event = 0; current_event < num_of_events; current_event++) {
    if (_loadbar && current_event % 10000 == 0) DataHandeler::loadbar(current_event, num_of_events);
    total += DataHandeler::Run(current_event);
  }
  _chain->Reset();
  return total;
}

int DataHandeler::Run(int current_event) {
  _chain->GetEntry(current_event);
  auto check = std::make_unique<Cuts>(_data);
  if (!check->isElecctron()) return 0;

  _hists->EC_inout(_data->ec_ei(0), _data->ec_eo(0));
  _hists->EC_fill(_data->etot(0), _data->p(0));
  _hists->TM_Fill(_data->p(0), physics::theta_calc(_data->cz(0)));
  double theta_cc = TMath::ACos(TMath::Abs(_data->p(0) * _data->cz(0)) / TMath::Abs(_data->p(0))) / D2R;
  _hists->CC_fill(_data->cc_sect(0), (_data->cc_segm(0) % 1000) / 10, _data->cc_segm(0) / 1000 - 1, _data->nphe(0),
                  theta_cc);

  auto event = std::make_unique<Reaction>(_data);

  _hists->Fill_E_Prime_fid(event->e_mu_prime());
  _hists->Fill_E_Prime(event->e_mu_prime());

  if (getenv("CUTS") != nullptr && atoi(getenv("CUTS")) == true) {
    CUTS = check->isStrictElecctron();
  } else {
    CUTS = check->isElecctron();
  }

  if (CUTS) {
    _hists->Fill_Beam_Position(_data->dc_vx(0), _data->dc_vy(0), _data->dc_vz(0));

    auto dt = std::make_unique<Delta_T>(_data->sc_t(0), _data->sc_r(0));

    dt->delta_t_hists(_hists, _data);
    std::vector<double> dt_proton = dt->delta_t_array(MASS_P, _data);
    std::vector<double> dt_pi = dt->delta_t_array(MASS_PIP, _data);

    float theta = physics::theta_calc(_data->cz(0));
    float phi = physics::phi_calc(_data->cx(0), _data->cy(0));
    int sector = _data->dc_sect(0);
    _hists->Fill_electron_fid(theta, phi, sector);

    // auto photon_flux = std::make_unique<PhotonFlux>(event->e_mu(), event->e_mu_prime());
    //_hists->Photon_flux_Fill(photon_flux->GetVirtualPhotonFlux());

    _hists->WvsQ2_Fill(event->W(), event->Q2(), _data->ec_sect(0));

    for (int part_num = 1; part_num < _data->gpart(); part_num++) {
      theta = physics::theta_calc(_data->cz(part_num));
      phi = physics::phi_calc(_data->cx(part_num), _data->cy(part_num));
      sector = _data->dc_sect(part_num);

      _hists->delta_t_sec_pad(_data->p(part_num), _data->q(part_num), dt->Get_dt_P(), dt->Get_dt_Pi(), dt->Get_dt_E(),
                              _data->sc_sect(part_num), _data->sc_pd(part_num));

      _hists->Fill_Target_Vertex(_data->vx(part_num), _data->vy(part_num), _data->vz(part_num));
      _hists->MomVsBeta_Fill(_data->p(part_num), _data->b(part_num));

      if (_data->q(part_num) == POSITIVE) {
        _hists->MomVsBeta_Fill_pos(_data->p(part_num), _data->b(part_num));
        if (check->dt_Pip_cut(dt_pi.at(part_num), _data->p(part_num)) &&
            check->dt_P_cut(dt_proton.at(part_num), _data->p(part_num)))
          _hists->Fill_proton_Pi_ID_P(_data->p(part_num), _data->b(part_num));

        if (check->dt_Pip_cut(dt_pi.at(part_num), _data->p(part_num))) {
          event->SetPip(part_num);
          _hists->Fill_hadron_fid(theta, phi, sector, PIP);
          _hists->Fill_pion_WQ2(event->W(), event->Q2());
          _hists->Fill_Pi_ID_P(_data->p(part_num), _data->b(part_num));
        } else if (check->dt_P_cut(dt_proton.at(part_num), _data->p(part_num))) {
          event->SetProton(part_num);
          _hists->Fill_hadron_fid(theta, phi, sector, PROTON);
          _hists->Fill_proton_WQ2(event->W(), event->Q2());
          _hists->Fill_proton_ID_P(_data->p(part_num), _data->b(part_num));
        } else
          event->SetOther(part_num);

      } else if (_data->q(part_num) == NEGATIVE) {
        _hists->MomVsBeta_Fill_neg(_data->p(part_num), _data->b(part_num));
        if (check->dt_Pip_cut(dt_pi.at(part_num), _data->p(part_num))) {
          _hists->Fill_hadron_fid(theta, phi, sector, PIM);
          event->SetPim(part_num);
        } else
          event->SetOther(part_num);
      } else if (_data->q(part_num) == 0) {
        _hists->MomVsBeta_Fill_neutral(_data->p(part_num), _data->b(part_num));
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
      _hists->Fill_channel_WQ2(event->W(), event->Q2(), _data->ec_sect(0), event->e_mu_prime(), event->MM(),
                               event->MM2());
      _hists->Fill_Missing_Mass_strict(event->MM(), event->MM2());
      _hists->EC_cut_fill(_data->etot(0), _data->p(0));
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

  _chain = new TChain("h10");
  for (auto& f : _input_files) _chain->Add(f.c_str());
  _data = std::make_shared<Branches>(_chain, true);
  int num_of_events = (int)_chain->GetEntries();
  int current_event = 0;

  for (current_event = 0; current_event < num_of_events; current_event++) {
    if (_loadbar && current_event % 10000 == 0) DataHandeler::loadbar(current_event, num_of_events);
    total += DataHandeler::Run(current_event);
    total += mcHandeler::Run(current_event);
  }
  _chain->Reset();
  return total;
}

int mcHandeler::Run(int current_event) {
  _chain->GetEntry(current_event);
  auto check = std::make_unique<Cuts>(_data);

  auto mc_event = std::make_unique<Reaction>(_data, true);
  _mc_hists->Fill_P(_data);
  _mc_hists->Fill_WQ2_MC(mc_event->W(), mc_event->Q2());

  for (int part_num = 1; part_num < _data->gpart(); part_num++) {
    if (_data->pidpart(part_num) == PIP) mc_event->SetPip(part_num);
  }
  if (mc_event->SinglePip()) _mc_hists->Fill_Missing_Mass(mc_event->MM(), mc_event->MM2());

  return 1;
}
