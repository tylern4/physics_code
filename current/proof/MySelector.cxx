#define MySelector_cxx
// The class definition in MySelector.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("MySelector.C")
// root> T->Process("MySelector.C","some options")
// root> T->Process("MySelector.C+")
//

#include <TH2.h>
#include <TStyle.h>
#include "MySelector.h"
void MySelector::Begin(TTree* /*tree*/) {
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();
}

void MySelector::SlaveBegin(TTree* /*tree*/) {
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

  fWq2 = new TH2D("fWq2", "W vs Q^{2}", 500, 0.0, 2.0, 500, 0.0, 5.0);
  fWq2->SetDirectory(0);
  fWq2->GetXaxis()->SetTitle("W");
  fWq2->GetYaxis()->SetTitle("Q^{2}");
  fOutput->Add(fWq2);

  fW = new TH1D("fW", "W", 500, 0.0, 2.0);
  fW->SetDirectory(0);
  fW->GetXaxis()->SetTitle("W");
  fOutput->Add(fW);
}

Bool_t MySelector::Process(Long64_t entry) {
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // When processing keyed objects with PROOF, the object is already loaded
  // and is available via the fObject pointer.
  //
  // This function should contain the \"body\" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.

  fReader.SetEntry(entry);

  bool electron_cuts = true;
  // electron cuts
  electron_cuts &= (id[0] == 11);             // First particle is electron
  electron_cuts &= (stat[0] > 0);             // First Particle hit stat
  electron_cuts &= ((int)q[0] == -1);         // First particle is negative Q
  electron_cuts &= (sc[0] > 0);               // First Particle hit sc
  electron_cuts &= (dc[0] > 0);               // ``` ``` ``` d
  electron_cuts &= (ec[0] > 0);               // ``` ``` ``` ec
  electron_cuts &= (dc_stat[dc[0] - 1] > 0);  //??
  electron_cuts &= (etot[ec[0] - 1] / p[0]) < 0.4;
  electron_cuts &= (etot[ec[0] - 1] / p[0]) > 0.2;

  if (electron_cuts) {
    TLorentzVector e_mu_prime;
    e_mu_prime.SetXYZM(p[0] * cx[0], p[0] * cy[0], p[0] * cz[0], 0.000511);
    fWq2->Fill(W_calc(e_mu_prime), Q2_calc(e_mu_prime));
    fW->Fill(W_calc(e_mu_prime));
    return kTRUE;
  } else {
    return kFALSE;
  }
}

void MySelector::SlaveTerminate() {
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
}

void MySelector::Terminate() {
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  TCanvas* c1 = new TCanvas("c1", "WvsQ2 Canvas", 1600, 900);
  c1->Divide(2);
  TH2D* hf_wq2 = dynamic_cast<TH2D*>(fOutput->FindObject("fWq2"));
  c1->cd(1);
  hf_wq2->Draw("colz");
  c1->Update();

  TH1D* hf_w = dynamic_cast<TH1D*>(fOutput->FindObject("fW"));
  c1->cd(2);
  hf_w->Draw();
  c1->Update();
}
