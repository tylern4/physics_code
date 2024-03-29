#define sim_cxx
// The class definition in sim.h has been generated automatically
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
// root> T->Process("sim.C")
// root> T->Process("sim.C","some options")
// root> T->Process("sim.C+")
//

#include <TH2.h>
#include <TStyle.h>
#include "sim.hh"

void sim::Begin(TTree* /*tree*/) {
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();
}

void sim::SlaveBegin(TTree* /*tree*/) {
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();

  for (int i = 0; i < 7; i++) {
    for (int j = 0; j < 4; j++) {
      mom_sec[i][j] = new TH1D(Form("P%c_sec%d", xyz[j], i), Form("#DeltaP_{%c} Sector %d", xyz[j], i), 500, -1, 1);
      fOutput->Add(mom_sec[i][j]);
    }
  }

  fWq2 = new TH2D("fWq2", "W vs Q^{2}", 500, 0.0, 3.0, 500, 0.0, 5.0);
  fWq2->SetDirectory(0);
  fWq2->GetXaxis()->SetTitle("W");
  fWq2->GetYaxis()->SetTitle("Q^{2}");
  fOutput->Add(fWq2);

  fW = new TH1D("fW", "W", 500, 0.0, 3.0);
  fW->SetDirectory(0);
  fW->GetXaxis()->SetTitle("W");
  fOutput->Add(fW);
}

Bool_t sim::Process(Long64_t entry) {
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
  TLorentzVector e_mu_prime;
  bool electron_cuts = true;
  // electron cuts
  // electron_cuts &= (id[0] == 11);             // First particle is electron
  electron_cuts &= (stat[0] > 0);             // First Particle hit stat
  electron_cuts &= ((int)q[0] == -1);         // First particle is negative Q
  electron_cuts &= (sc[0] > 0);               // First Particle hit sc
  electron_cuts &= (dc[0] > 0);               // ``` ``` ``` d
  electron_cuts &= (ec[0] > 0);               // ``` ``` ``` ec
  electron_cuts &= (dc_stat[dc[0] - 1] > 0);  //??

  e_mu_prime.SetXYZM(p[0] * cx[0], p[0] * cy[0], p[0] * cz[0], 0.000511);

  if (electron_cuts) {
    fWq2->Fill(W_calc(e_mu_prime), Q2_calc(e_mu_prime));
    fW->Fill(W_calc(e_mu_prime));

    mom_sec[0][0]->Fill((pxpart[0] - p[0] * cx[0]) / pxpart[0]);
    mom_sec[0][1]->Fill((pypart[0] - p[0] * cy[0]) / pypart[0]);
    mom_sec[0][2]->Fill((pzpart[0] - p[0] * cz[0]) / pzpart[0]);
    mom_sec[0][3]->Fill((TMath::Sqrt(pxpart[0] * pxpart[0] + pypart[0] * pypart[0] + pzpart[0] * pzpart[0]) - p[0]) /
                        TMath::Sqrt(pxpart[0] * pxpart[0] + pypart[0] * pypart[0] + pzpart[0] * pzpart[0]));

    int sec = ec_sect[0];
    if (sec < 6) {
      mom_sec[sec + 1][0]->Fill((pxpart[0] - p[0] * cx[0]) / pxpart[0]);
      mom_sec[sec + 1][1]->Fill((pypart[0] - p[0] * cy[0]) / pypart[0]);
      mom_sec[sec + 1][2]->Fill((pzpart[0] - p[0] * cz[0]) / pzpart[0]);
      mom_sec[sec + 1][3]->Fill(
          (TMath::Sqrt(pxpart[0] * pxpart[0] + pypart[0] * pypart[0] + pzpart[0] * pzpart[0]) - p[0]) /
          TMath::Sqrt(pxpart[0] * pxpart[0] + pypart[0] * pypart[0] + pzpart[0] * pzpart[0]));
    }
  }

  return kTRUE;
}

void sim::SlaveTerminate() {
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
}

void sim::Terminate() {
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  TFile* outfile = new TFile("out.root", "RECREATE");
  outfile->cd();

  for (int i = 0; i < 7; i++) {
    for (int j = 0; j < 4; j++) {
      TH1D* hf = dynamic_cast<TH1D*>(fOutput->FindObject(Form("P%c_sec%d", xyz[j], i)));
      hf->Write();
    }
  }

  TCanvas* c1 = new TCanvas("wq2", "WvsQ2 Canvas", 1600, 900);
  c1->Divide(2);
  TH2D* hf_wq2 = dynamic_cast<TH2D*>(fOutput->FindObject("fWq2"));
  c1->cd(1);
  hf_wq2->Draw("colz");

  TH1D* hf_w = dynamic_cast<TH1D*>(fOutput->FindObject("fW"));
  c1->cd(2);
  hf_w->Draw();

  hf_wq2->Write();
  hf_w->Write();
  c1->Write();
}
