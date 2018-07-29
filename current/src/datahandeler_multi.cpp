/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#include "datahandeler_multi.hpp"

Datahandeler_multi::Datahandeler_multi(std::vector<std::string> fin, std::string RootFile_output) {
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
}

Datahandeler_multi::~Datahandeler_multi() {}

void Datahandeler_multi::loadbar(long x, long n) {
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

void Datahandeler_multi::run() {
  int size = input_files.size();
  int i = 0;
  std::thread *fh_thread[size];

#ifdef __PID__
  std::cerr << CYAN << "Using builtin PID" << DEF << std::endl;
#endif

  for (i = 0; i < size; i++) {
    loadbar(i, size - 1);
#ifndef __THREAD__
    load_files(input_files.at(i));
#endif

#ifdef __THREAD__
    try {
      fh_thread[i] = new std::thread(std::mem_fn(&Datahandeler_multi::load_files), this, input_files.at(i));
      fh_thread[i]->join();
    } catch (const std::exception &e) {
      std::cerr << RED << "Error:\t" << e.what() << std::endl;
      std::cerr << CYAN << "Bad File: \t" << input_files.at(i) << DEF << std::endl;
    }
#endif
  }

  int tot_num_events = All_events.size();
  int tot, p;
  for (tot = 0; tot < tot_num_events; tot++) {
    int num_parts = All_events[tot].size();
    TLorentzVector reaction = TLorentzVector(0.0, 0.0, 0.0, MASS_P);
    for (p = 0; p < num_parts; p++) {
      TLorentzVector particle = All_events[tot][p];
      if (tot == 0) {  // electron
        TLorentzVector gamma_mu = (*e_mu - particle);
        reaction += gamma_mu;
      } else {  // particles
        reaction -= particle;
      }
    }
    mm->Fill(reaction.Mag2());
  }
  mm->Write();
}

void Datahandeler_multi::load_files(std::string fin) {
  // Load chain from branch h10
  std::vector<TLorentzVector> Events;
  bool cuts, electron_cuts;
  int num_of_events;
  int total_events;
  int num_of_pis;
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
    check->Set_BeamPosition((double)dc_vx[dc[0] - 1], (double)dc_vy[dc[0] - 1], (double)dc_vz[dc[0] - 1]);

    if (check->isStrictElecctron()) {
      Events.resize(gpart);
      TLorentzVector e_mu_prime = physics::fourVec(p[0], cx[0], cy[0], cz[0], MASS_E);
      Events[0] = e_mu_prime;
      TLorentzVector gamma_mu = (*e_mu - e_mu_prime);

      Delta_T *dt = new Delta_T(sc_t[sc[0] - 1], sc_r[sc[0] - 1]);
      std::vector<double> dt_proton = dt->delta_t_array(MASS_P, gpart);
      std::vector<double> dt_pi = dt->delta_t_array(MASS_PIP, gpart);
      int num_pip = 0;
      int num_p = 0;
      for (int part_num = 1; part_num < gpart; part_num++) {
        if (q[part_num] == POSITIVE) {
          if (check->dt_P_cut(dt_proton.at(part_num), p[part_num])) {
            PID = PROTON;
            num_p++;
          }
          if (check->dt_Pip_cut(dt_proton.at(part_num), p[part_num])) {
            PID = PIP;
            num_pip++;
          }
        } else if (q[part_num] == NEGATIVE) {
          if (check->dt_Pip_cut(dt_proton.at(part_num), p[part_num])) PID = PIM;
        }
        TLorentzVector Particle = physics::fourVec(p[part_num], cx[part_num], cy[part_num], cz[part_num], PID);
        Events[part_num] = Particle;
      }
      if (num_pip == 1) All_events.emplace_back(Events);
      Events.clear();
    }
    delete check;
  }
  chain->Reset();  // delete Tree object
}
