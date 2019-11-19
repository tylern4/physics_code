/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/
#include "datahandeler.hpp"

DataHandeler::DataHandeler() = default;

DataHandeler::DataHandeler(const std::vector<std::string>& fin, const std::shared_ptr<Histogram>& hists)
    : _input_files(fin), _hists(hists) {
  _chain = std::make_shared<TChain>("h10");
  for (auto& f : _input_files) _chain->Add(f.c_str());
  _data = std::make_shared<Branches>(_chain);
}

DataHandeler::DataHandeler(const std::vector<std::string>& fin, const std::shared_ptr<Histogram>& hists,
                           const std::shared_ptr<MomCorr>& mom_corr)
    : _input_files(fin), _hists(hists), _mom_corr(mom_corr) {
  _chain = std::make_shared<TChain>("h10");
  for (auto& f : _input_files) _chain->Add(f.c_str());
  _data = std::make_shared<Branches>(_chain);
}

DataHandeler::~DataHandeler() = default;

void DataHandeler::loadbar(long x, long n) {
  int w = 50;
  if (x >= n) return;

  double ratio = x / (double)n;
  int c = ratio * w;

  std::cerr << BLUE << " [";
  for (int x = 0; x < c; x++) std::cerr << GREEN << "=" << DEF;
  std::cerr << GREEN << ">" << DEF;
  for (int x = c; x < w; x++) std::cerr << " ";
  std::cerr << BLUE << (int)(ratio * 100) << "%]\r" << DEF << std::flush;
}

void DataHandeler::setLoadBar(bool load) { _loadbar = load; }

int DataHandeler::Run() {
  if (_hists == nullptr) return 0;
  size_t num_of_events = (size_t)_chain->GetEntries();
  for (size_t current_event = 0; current_event < num_of_events; current_event++) {
    if (_loadbar && current_event % 10000 == 0) DataHandeler::loadbar(current_event, num_of_events);
    DataHandeler::RunEvent(current_event);
  }
  _chain->Reset();
  return num_of_events;
}

const void DataHandeler::RunEvent(size_t current_event) {
  _chain->GetEntry(current_event);
  auto check = std::make_unique<Cuts>(_data);
  if (!check->isElecctron()) return;

  _hists->EC_inout(_data->ec_ei(0), _data->ec_eo(0));
  _hists->EC_fill(_data->etot(0), _data->p(0));
  _hists->TM_Fill(_data->p(0), physics::theta_calc(_data->cz(0)));
  double theta_cc = TMath::ACos(TMath::Abs(_data->p(0) * _data->cz(0)) / TMath::Abs(_data->p(0))) / D2R;
  _hists->CC_fill(_data->cc_sect(0), (_data->cc_segm(0) % 1000) / 10, _data->cc_segm(0) / 1000 - 1, _data->nphe(0),
                  theta_cc);

  auto event = std::make_shared<Reaction>(_data, _mom_corr);

  _hists->Fill_E_Prime_fid(event->e_mu_prime());
  _hists->Fill_E_Prime(event->e_mu_prime());

  if (getenv("CUTS") != nullptr && atoi(getenv("CUTS")) == true) {
    CUTS = check->isStrictElecctron();
  } else {
    CUTS = check->isElecctron();
  }

  if (!CUTS) return;
  _hists->Fill_Beam_Position(_data->dc_vx(0), _data->dc_vy(0), _data->dc_vz(0));

  auto photon_flux = std::make_unique<PhotonFlux>(event->e_mu(), event->e_mu_prime());
  _hists->Photon_flux_Fill(photon_flux->GetVirtualPhotonFlux());
  _hists->WvsQ2_Fill(event->W(), event->Q2(), _data->ec_sect(0));

  auto dt = std::make_unique<Delta_T>(_data);
  dt->delta_t_hists(_hists);

  float theta = physics::theta_calc(_data->cz(0));
  float phi = physics::phi_calc(_data->cx(0), _data->cy(0));
  int sector = _data->dc_sect(0);

  _hists->Fill_electron_fid(theta, phi, sector);

  for (int part_num = 1; part_num < _data->gpart(); part_num++) {
    theta = physics::theta_calc(_data->cz(part_num));
    phi = physics::phi_calc(_data->cx(part_num), _data->cy(part_num));
    sector = _data->dc_sect(part_num);

    _hists->delta_t_sec_pad(_data->p(part_num), _data->q(part_num), dt->Get_dt_P(part_num), dt->Get_dt_Pi(part_num),
                            dt->Get_dt_E(part_num), _data->sc_sect(part_num), _data->sc_pd(part_num));

    _hists->Fill_Target_Vertex(_data->vx(part_num), _data->vy(part_num), _data->vz(part_num));
    _hists->MomVsBeta_Fill(_data->p(part_num), _data->b(part_num));

    if (_data->q(part_num) == POSITIVE) {
      _hists->MomVsBeta_Fill_pos(_data->p(part_num), _data->b(part_num));
    } else if (_data->q(part_num) == NEGATIVE) {
      _hists->MomVsBeta_Fill_neg(_data->p(part_num), _data->b(part_num));
    } else if (_data->q(part_num) == 0) {
      if (_data->id(part_num) == NEUTRON) {
        event->SetNeutron(part_num);
        _hists->MomVsBeta_Fill_neutral(_data->p(part_num), _data->b(part_num));
        _hists->Fill_neutron_fid(_data->cc_c2(part_num), _data->cc_r(part_num), _data->cc_sect(part_num));
      } else
        event->SetOther(part_num);
      continue;
    }

    if (check->Pip(part_num) && check->Prot(part_num))
      _hists->Fill_proton_Pi_ID_P(_data->p(part_num), _data->b(part_num));

    if (check->Pip(part_num)) {
      event->SetPip(part_num);
      _hists->Fill_hadron_fid(theta, phi, sector, PIP);
      _hists->Fill_pion_WQ2(event->W(), event->Q2());
      _hists->Fill_Pi_ID_P(_data->p(part_num), _data->b(part_num));
    } else if (check->Prot(part_num)) {
      event->SetProton(part_num);
      _hists->Fill_hadron_fid(theta, phi, sector, PROTON);
      _hists->Fill_proton_WQ2(event->W(), event->Q2());
      _hists->Fill_proton_ID_P(_data->p(part_num), _data->b(part_num));
    } else if (check->Pim(part_num)) {
      _hists->Fill_hadron_fid(theta, phi, sector, PIM);
      event->SetPim(part_num);
    } else
      event->SetOther(part_num);
  }

  if (event->channel()) _hists->EC_cut_fill(_data->etot(0), _data->p(0));
  _hists->FillEvent(event);

  return;
}

mcHandeler::mcHandeler(const std::vector<std::string>& fin, const std::shared_ptr<mcHistogram>& hists) {
  _input_files = fin;
  _hists = hists;
  _mc_hists = hists;
  _chain = std::make_shared<TChain>("h10");
  for (auto& f : _input_files) _chain->Add(f.c_str());
  _data = std::make_shared<Branches>(_chain, true);
}

int mcHandeler::Run() {
  size_t num_of_events = (size_t)_chain->GetEntries();
  for (size_t current_event = 0; current_event < num_of_events; current_event++) {
    if (_loadbar && current_event % 10000 == 0) DataHandeler::loadbar(current_event, num_of_events);
    DataHandeler::RunEvent(current_event);
    mcHandeler::RunEvent(current_event);
  }
  _chain->Reset();
  return num_of_events;
}

const void mcHandeler::RunEvent(int current_event) {
  _chain->GetEntry(current_event);
  auto check = std::make_unique<Cuts>(_data);

  auto mc_event = std::make_shared<MCReaction>(_data);
  _mc_hists->Fill_P(_data);
  _mc_hists->Fill_WQ2_MC(mc_event);
  _hists->Fill_ND(mc_event);
  _mc_hists->Fill(mc_event);

  for (int part_num = 1; part_num < _data->gpart(); part_num++) {
    if (_data->pidpart(part_num) == PIP) mc_event->SetPip(part_num);
  }
  if (mc_event->SinglePip()) _mc_hists->Fill_Missing_Mass(mc_event);

  return;
}
