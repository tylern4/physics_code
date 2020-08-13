/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#include "yeilds.hpp"

Yeilds::Yeilds() {}
Yeilds::Yeilds(std::string output_file_name) { csv_output.open(output_file_name); }
Yeilds::Yeilds(std::string output_file_name, bool isRoot = true) {
  Rootout = std::make_shared<TFile>(output_file_name.c_str(), "RECREATE");
  ntuple =
      std::make_shared<TNtuple>("ntuple", "", "type:W:Q2:MM:MM2:theta_e:theta_star:phi_star:theta_lab:phi_lab:sector");
}
Yeilds::~Yeilds() {
  if (ntuple) ntuple->Write();
}

void Yeilds::OpenFile(std::string output_file_name) { csv_output.open(output_file_name); }

void Yeilds::WriteHeader() {
  csv_output << "electron_sector,w,q2,theta,phi,mm2,e_p,e_cx,e_cy,e_cz,pip_p,pip_cx,pip_cy,pip_cz,helicty,type"
             << std::endl;
}
