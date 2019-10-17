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

std::string mom_correction_csv(const std::vector<std::string>& fins) {
  std::string mom_correction;
  auto chain = std::make_shared<TChain>("h10");
  for (auto& fin : fins) chain->Add(fin.c_str());
  auto data = std::make_shared<Branches>(chain);

  size_t num_of_events = chain->GetEntries();

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain->GetEntry(current_event);
    auto cuts = std::make_unique<Cuts>(data);

    if (!cuts->isElecctron()) continue;
    if (!cuts->Beam_cut()) continue;

    auto event = std::make_shared<Reaction>(data);
    int prot_num = -1;

    for (int part_num = 1; part_num < data->gpart(); part_num++) {
      if (cuts->Prot(part_num)) {
        prot_num = part_num;
        event->SetProton(part_num);
      } else
        event->SetOther(part_num);
    }
    if (event->SingleP() && physics::theta_calc(data->pz(prot_num)) > 30 && event->W() > 0.5 && event->W() < 1.1)
      mom_correction += std::to_string(data->px(0)) + "," + std::to_string(data->py(0)) + "," +
                        std::to_string(data->pz(0)) + "," + std::to_string(data->px(prot_num)) + "," +
                        std::to_string(data->py(prot_num)) + "," + std::to_string(data->pz(prot_num)) + "," +
                        std::to_string(event->W()) + "," + std::to_string(event->Q2()) + "," +
                        std::to_string(event->sector()) + "\n";
  }

  return mom_correction;
}
#endif
