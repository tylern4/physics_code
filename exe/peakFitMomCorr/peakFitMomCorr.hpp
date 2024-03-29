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

std::string mom_correction_csv(const std::vector<std::string>& fins, bool MC, size_t thread) {
  std::string mom_correction;
  auto chain = std::make_shared<TChain>("h10");
  for (auto& fin : fins) chain->Add(fin.c_str());
  auto data = std::make_shared<Branches>(chain, MC);

  size_t num_of_events = chain->GetEntries();

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    if (thread == 0 && current_event % 1000 == 0)
      std::cerr << "\t" << 100 * current_event / num_of_events << "\r" << std::flush;
    chain->GetEntry(current_event);
    auto cuts = std::make_unique<Cuts>(data);

    auto e_mc = LorentzVector(data->pxpart(0), data->pypart(0), data->pzpart(0), mass_map[ELECTRON]);
    mom_correction += std::to_string(e_mc.P()) + "," + std::to_string(e_mc.Theta()) + "," + std::to_string(e_mc.Phi());
    if (!cuts->isElectron()) {
      mom_correction += ",,,,,None\n";
      continue;
    }
    if (!cuts->Beam_cut()) {
      mom_correction += ",,,,,None\n";
      continue;
    }

    auto event = std::make_shared<Reaction>(data, nullptr);
    float theta = physics::theta_calc(data->cz(0));
    float phi = physics::phi_calc(data->cx(0), data->cy(0));
    short sector = data->dc_sect(0);

    for (int part_num = 1; part_num < data->gpart(); part_num++) {
      if (cuts->Pip(part_num))
        event->SetPip(part_num);
      else if (cuts->Prot(part_num))
        event->SetProton(part_num);
      else if (cuts->Pim(part_num))
        event->SetPim(part_num);
      else
        event->SetOther(part_num);
    }

    std::string type = "unknown";
    if (event->channel())
      type = "Npip";
    else if (event->elastic())
      type = "elastic";
    else
      type = "";

    if (type != "")
      mom_correction += "," + std::to_string(data->p(0)) + "," + std::to_string(theta * DEG2RAD) + "," +
                        std::to_string(phi * DEG2RAD) + "," + std::to_string(event->sector()) + "," + type + "\n";
    else {
      mom_correction += ",,,,,None\n";
    }
  }

  return mom_correction;
}
#endif
