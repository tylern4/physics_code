/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef MOM_CORRECTIONS_H
#define MOM_CORRECTIONS_H
#include <string>
#include "branches.hpp"
#include "constants.hpp"
#include "cuts.hpp"
#include "syncfile.hpp"

size_t mom_correction_csv(const std::vector<std::string>& fins, const std::shared_ptr<SyncFile>& sync,
                          const std::shared_ptr<MomCorr>& mom_corr) {
  std::string mom_correction;
  auto chain = std::make_shared<TChain>("h10");
  for (auto& fin : fins) chain->Add(fin.c_str());
  auto data = std::make_shared<Branches>(chain);

  size_t num_of_events = chain->GetEntries();

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain->GetEntry(current_event);
    auto cuts = std::make_unique<e1d_Cuts>(data);

    if (!cuts->check_banks()) continue;
    if (!cuts->isElectron()) continue;
    if (!cuts->Beam_cut()) continue;

    auto event = std::make_shared<Reaction>(data, mom_corr);
    int positive_num = -1;

    for (int part_num = 1; part_num < data->gpart(); part_num++) {
      if (cuts->Pip(part_num)) {
        positive_num = part_num;
        event->SetPip(part_num);
      } else if (cuts->Prot(part_num)) {
        positive_num = part_num;
        event->SetProton(part_num);
      } else if (cuts->Pim(part_num)) {
        event->SetPim(part_num);
      } else
        event->SetOther(part_num);
    }

    if (event->elastic() || true) {
      std::string type = "none";
      if (event->elastic()) type = "elastic";
      if (event->channel()) type = "channel";

      auto mom_correction = std::to_string(data->p(0)) + "," + std::to_string(physics::theta_calc_rad(data->cz(0))) +
                            "," + std::to_string(physics::phi_calc_rad(data->cx(0), data->cy(0))) + "," +
                            std::to_string(data->p(positive_num)) + "," +
                            std::to_string(physics::theta_calc_rad(data->cz(positive_num))) + "," +
                            std::to_string(physics::phi_calc_rad(data->cx(positive_num), data->cy(positive_num))) +
                            "," + std::to_string(event->W()) + "," + std::to_string(event->Q2()) + "," +
                            std::to_string(event->sector()) + "," + type + "\n";
      sync->write(mom_correction);
    }
  }
  return num_of_events;
}
#endif
