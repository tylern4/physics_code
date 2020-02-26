#include <iostream>
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "THn.h"
#include "THnSparse.h"
#include "TLeaf.h"
#include "TStyle.h"
#include "TTree.h"

void branch_add_pfermi(std::string root_file, std::string weights_file) {
  TFile *file2 = TFile::Open(weights_file.c_str(), "READ");
  TFile *file1 = TFile::Open(root_file.c_str(), "UPDATE");

  TTree *T2 = (TTree *)file2->Get("h10");

  TBranch *br_sig = T2->GetBranch("sigma");
  TBranch *br_mom = T2->GetBranch("p_el_test");
  TBranch *br_px_f = T2->GetBranch("px_fermi");
  TBranch *br_py_f = T2->GetBranch("py_fermi");
  TBranch *br_pz_f = T2->GetBranch("pz_fermi");

  // TBranch *newBranch = t3->Branch("new_v", &new_v, "new_v/F");

  TTree *T1 = (TTree *)file1->Get("h10");
  Float_t sigma_tot, p_el_new2, px_ferm, py_ferm, pz_ferm;

  //  T1->Branch("sigma1",&sigma_total,"sigma1/F");
  //  T1->Branch("p_el_test1",&p_el_test2,"p_el_test1/F");

  TBranch *sigma_total = T1->Branch("sigma_total", &sigma_tot, "sigma_total/F");
  TBranch *p_el_new = T1->Branch("p_el_new", &p_el_new2, "p_el_new/F");
  TBranch *px_Fermi = T1->Branch("px_Fermi", &px_ferm, "px_Fermi/F");
  TBranch *py_Fermi = T1->Branch("py_Fermi", &py_ferm, "py_Fermi/F");
  TBranch *pz_Fermi = T1->Branch("pz_Fermi", &pz_ferm, "pz_Fermi/F");

  for (int i = 0; i < br_sig->GetEntries(); i++) {
    br_sig->GetEntry(i);
    br_mom->GetEntry(i);
    br_px_f->GetEntry(i);
    br_py_f->GetEntry(i);
    br_pz_f->GetEntry(i);

    sigma_tot = br_sig->GetLeaf("sigma")->GetValue();
    p_el_new2 = br_mom->GetLeaf("p_el_test")->GetValue();

    px_ferm = br_px_f->GetLeaf("px_fermi")->GetValue();
    py_ferm = br_py_f->GetLeaf("py_fermi")->GetValue();
    pz_ferm = br_pz_f->GetLeaf("pz_fermi")->GetValue();

    T1->SetDirectory(0);
    sigma_total->Fill();
    p_el_new->Fill();
    px_Fermi->Fill();
    py_Fermi->Fill();
    pz_Fermi->Fill();
  }
  T1->Write("", TObject::kOverwrite);
}

#if not defined(__CLING__)
int main(int argc, char const *argv[]) {
  if (argc < 2) {
    std::cerr << "Not enough arguments" << std::endl;
    std::cerr << "To Use:\t" << argv[0] << " data.root mc.root" << std::endl;
    exit(1);
  }

  branch_add_pfermi(argv[1], argv[2]);

  return 0;
}
#endif