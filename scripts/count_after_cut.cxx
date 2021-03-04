#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "TCanvas.h"
#include "TMath.h"
//#include "vector"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TTree.h"

using namespace std;
void count_after_cut(std::string filename_to_open) {
  // Set constant variable for mass
  ifstream infile(filename_to_open);
  string s, str_runnum, str_filenum;
  TFile *rootfile;
  TTree *tree;

  if (!infile) {
    cout << "Can't open infile " << endl;
    exit(1);
  }

  ofstream outfile("count_event", ios::out);

  while (infile >> s) {
    if (s == "end") {
      cout << "end of list" << endl;
      break;
    }
    // cout << s << endl;
    rootfile = new TFile(s.c_str());
    tree = (TTree *)rootfile->Get("h10");
    // Set variables for branches
    // Trip_flag: identify if events are good or bad
    Float_t q_l;
    Float_t t_l;
    Float_t tr_time;
    Float_t rf_time1;
    // Int_t latch1;
    // Int_t hlsc;
    Int_t intt;
    Int_t gpart;
    Int_t id[20];    //[gpart]
    Int_t stat[20];  //[gpart]
    Int_t dc[20];    //[gpart]
    Int_t cc[20];    //[gpart]
    Int_t sc[20];    //[gpart]
    Int_t ec[20];    //[gpart]
    Int_t lec[20];   //[gpart]
    Float_t p[20];   //[gpart]
    Int_t q[20];     //[gpart]
    Float_t b[20];   //[gpart]
    Float_t cx[20];  //[gpart]
    Float_t cy[20];  //[gpart]
    Float_t cz[20];  //[gpart]
    Float_t vx[20];  //[gpart]
    Float_t vy[20];  //[gpart]
    Float_t vz[20];  //[gpart]

    tree->SetBranchAddress("q_l", &q_l);
    tree->SetBranchAddress("t_l", &t_l);
    tree->SetBranchAddress("tr_time", &tr_time);
    tree->SetBranchAddress("rf_time1", &rf_time1);
    // tree->SetBranchAddress("latch1", &latch1);
    // tree->SetBranchAddress("hlsc", &hlsc);
    tree->SetBranchAddress("intt", &intt);
    tree->SetBranchAddress("gpart", &gpart);
    tree->SetBranchAddress("id", id);
    tree->SetBranchAddress("stat", stat);
    tree->SetBranchAddress("dc", dc);
    tree->SetBranchAddress("cc", cc);
    tree->SetBranchAddress("sc", sc);
    tree->SetBranchAddress("ec", ec);
    tree->SetBranchAddress("lec", lec);
    tree->SetBranchAddress("p", p);
    tree->SetBranchAddress("q", q);
    tree->SetBranchAddress("b", b);
    tree->SetBranchAddress("cx", cx);
    tree->SetBranchAddress("cy", cy);
    tree->SetBranchAddress("cz", cz);
    tree->SetBranchAddress("vx", vx);
    tree->SetBranchAddress("vy", vy);
    tree->SetBranchAddress("vz", vz);

    // Define variables
    float totalQ = 0;
    float qcurr;
    float qprev;
    float deltaq;
    float q_temp;
    TLorentzVector ni_4vec(0, 0, 0, 0);
    TVector3 ni_3vec(0, 0, 0);
    TLorentzVector n_4vec(0, 0, 0, 0);
    TVector3 n_3vec(0, 0, 0);
    TVector3 q_3vec(0, 0, 0);
    ni_4vec.SetVectM(ni_3vec, 0.9396);
    TLorentzVector ei_4vec(0, 0, 0, 0);
    TVector3 ei_3vec(0, 0, 2.039);
    ei_4vec.SetVectM(ei_3vec, 0.000511);
    TVector3 ef_3vec(0, 0, 0);
    TLorentzVector ef_4vec(0, 0, 0, 0);
    TVector3 pionf_3vec(0, 0, 0);
    TLorentzVector pionf_4vec(0, 0, 0, 0);
    TVector3 pf_3vec(0, 0, 0);
    TLorentzVector pf_4vec(0, 0, 0, 0);
    TVector3 protonf_3vec(0, 0, 0);
    TLorentzVector protonf_4vec(0, 0, 0, 0);
    TLorentzVector w_4vec(0, 0, 0, 0);
    TLorentzVector w_n_q_4vec(0, 0, 0, 0);
    TLorentzVector ZERO_4vec(0, 0, 0, 0);
    TLorentzVector pionf_CM_4vec(0, 0, 0, 0);
    // TLorentzVector cmq_4vec(0,0,0,0);
    // TLorentzVector cmPip_4vec;
    TVector3 b_3vec(0, 0, 0);
    TLorentzVector m_4vec(0, 0, 0, 0);

    Long64_t nentries = tree->GetEntries();

    int n_event = 0;
    for (Long64_t i = 0; i < nentries; i++) {
      ef_4vec = ZERO_4vec;
      pionf_4vec = ZERO_4vec;
      // Get the ith entry
      tree->GetEntry(i);
      // Get total farda cup charge
      q_temp = q_l;
      qcurr = q_temp;
      // cout<<"q_l="<<q_l<<endl;

      if (q_temp > 0.) {
        // cout<<"q_l"<<q_temp<<"qcurr"<<qcurr<<endl;
        if ((qcurr > qprev)) {
          deltaq = qcurr - qprev;
          totalQ += deltaq;
          // cout << setw(10) << "qcurr= " << qcurr << setw(10) << " qprev= " << qprev << setw(10) << " deltaq= " <<
          // deltaq
          //      << endl;
          // cout<<"q_l"<<q_temp<<"qcurr"<<qcurr<<"qprev"<<qprev<<"deltaq"<<deltaq<<endl;
        }
        qprev = qcurr;
      }
      // Cut through trip_flag

      for (int j = 0; j < gpart; j++) {
        // Cut delta time of good photons for SC && ST
        if (id[0] == 11) {
          ef_3vec.SetXYZ(p[0] * cx[0], p[0] * cy[0], p[0] * cz[0]);
          ef_4vec.SetVectM(ef_3vec, 0.000511);
          // cout<<"electron"<<endl;
        }

        if (id[j] == 211) {
          // if(dc[i]>0 && cc[i]>0 && sc[i]>0){
          pionf_3vec.SetXYZ(p[j] * cx[j], p[j] * cy[j], p[j] * cz[j]);
          pionf_4vec.SetVectM(pionf_3vec, 0.13957);
          // cout<<"pion"<<endl;
          // Pmom_pif_hist->Fill(p[i]);
        }
        if (ef_4vec != ZERO_4vec && pionf_4vec != ZERO_4vec) {
          n_event++;
        }
      }
    }
    tree->Delete();
    rootfile->Close();

    auto file_name = s.substr(s.find_last_of("/\\") + 1);
    str_runnum = file_name.substr(file_name.find("_r2") + 2, 5);
    str_filenum = file_name.substr(file_name.find("_r2") + 8, 2);

    // str_runnum = s.substr(s.length() - 24, 5);
    // str_filenum = s.substr(s.length() - 11, 2);

    cout << setw(10) << setiosflags(ios::left) << str_runnum << setw(10) << setiosflags(ios::left) << str_filenum
         << setw(10) << setiosflags(ios::left) << n_event << "totalQ=" << totalQ << endl;
    if (n_event != 0 && totalQ != 0) {
      outfile << str_runnum << "," << str_filenum << "," << n_event << "," << totalQ << "," << n_event / totalQ << endl;
    }
  }
  outfile.close();
}

#if !defined(__CLING__)
int main(int argc, char const *argv[]) {
  if (argc < 1) {
    std::cerr << "Not enough arguments" << std::endl;
    std::cerr << "To Use:\t" << argv[0] << " data.root" << std::endl;
    exit(1);
  }
  count_after_cut(argv[1]);
  return 0;
}
#endif