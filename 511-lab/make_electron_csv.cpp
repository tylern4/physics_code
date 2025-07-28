/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

// Only My Includes. All others in main.h
#include "511_lab.hpp"
#include "TStopwatch.h"
#include "branches.hpp"
#include "classes.hpp"
#include "constants.hpp"
#include "main.h"

void __make_electron_csv(std::string fin, std::string fout) {
  if (getenv("BEAM_E") != NULL) {
    BEAM_ENERGY = atof(getenv("BEAM_E"));
    std::cout << RED << "Beam energy set to: " << BEAM_ENERGY << DEF << std::endl;
  } else {
    BEAM_ENERGY = E1D_E0;
  }
  const char *progress = "-\\|/";
  int num_of_events;
  bool electron_cuts;
  double _p, _cx, _cy, _cz;

  auto chain = std::make_shared<TChain>("h10");
  chain->Add(fin.c_str());
  auto data = std::make_shared<Branches>(chain);

  // in main.h now
  // ofstream cut_outputs;
  csv_output.open(fout);
  num_of_events = (int)chain->GetEntries();

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain->GetEntry(current_event);
    if (current_event % 1000000 == 0)
      cout << "\t[ " << progress[((current_event / 1000000) % 4)] << " ]\t\t"
           << 100 * ((float)current_event / (float)num_of_events) << "\r\r" << flush;

    electron_cuts = true;
    // electron cuts
    electron_cuts &= (data->gpart() > 0);  // Number of good particles is gt 0
    electron_cuts &= (data->stat(0) > 0);  // First Particle hit stat
    electron_cuts &= (data->q(0) == -1);   // First particle is negative Q
    electron_cuts &= (data->sc(0) > 0);    // First Particle hit sc
    electron_cuts &= (data->dc(0) > 0);    // ``` ``` ``` dc
    electron_cuts &= (data->ec(0) > 0);    // ``` ``` ``` ec
    electron_cuts &= (data->dc_stat(0) > 0);

    if (!electron_cuts) continue;

    int n_prot = 0;
    int n_other = 0;
    for (int x = 1; x < data->gpart(); x++)
      if (data->id(x) == PROTON)
        n_prot++;
      else
        n_other++;

    if (n_prot == 1 && n_other == 0) {
      csv_output << data->p(0) << ",";
      csv_output << data->cx(0) << ",";
      csv_output << data->cy(0) << ",";
      csv_output << data->cz(0);
      csv_output << endl;
    }
  }

  csv_output.close();
}

int main(int argc, char **argv) {
  TStopwatch *Watch = new TStopwatch;
  Watch->Start();

  if (argc == 3) {
    std::string infilename = argv[1];
    std::string outfilename = argv[2];
    __make_electron_csv(infilename, outfilename);
  } else {
    std::cerr << RED << "Wrong inputs " << argv[0] << " input.root output.csv" << DEF << std::endl;
  }

  Watch->Stop();
  std::cout << RED << Watch->RealTime() << "sec" << DEF << std::endl;

  return 0;
}
