/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#include "yeilds.hpp"

Yeilds::Yeilds() {}

Yeilds::Yeilds(std::shared_ptr<SyncFile> multi_threaded_csv) : _multi_threaded_csv(multi_threaded_csv) {
  _multi_threaded_csv->write(Header());
}

Yeilds::Yeilds(std::string output_file_name, bool isRoot = true) {
  Rootout = std::make_shared<TFile>(output_file_name.c_str(), "RECREATE");
  ntuple = std::make_shared<TNtuple>("ntuple", "", "electron_sector:w:q2:theta:phi:mm2:helicty:type:hash");
}
Yeilds::~Yeilds() {
  if (ntuple) ntuple->Write();
}

std::string Yeilds::Header() { return "electron_sector,w,q2,theta,phi,mm2,helicty,type,hash"; }

void Yeilds::WriteData(const csv_data& toWrite) { _multi_threaded_csv->write(toWrite); };
