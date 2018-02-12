#define MyProof_cxx
// The class definition in MyProof.h has been generated automatically
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
// root> T->Process("MyProof.C")
// root> T->Process("MyProof.C","some options")
// root> T->Process("MyProof.C+")
//
#include "TROOT.h"
#include "TSystem.h"
#include "MyProof.h"
#include <TH2.h>
#include <TStyle.h>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"

static const double PI = TMath::Pi();
static const double D2R = PI / 180.0;
TH1D *mom;
TH2D *fiducial;
TH2D *fiducial_cut;
UInt_t fNumberOfEvents;
TDatime tBegin, tNow;
TVector3 e_mu_prime_3;
TLorentzVector e_mu_prime;

double fiducial_phi_lo(double theta_e, double theta_e_min, double k, double m) {
        return -37.14 * pow(sin((theta_e - theta_e_min) * 0.01745),
                            k + m / theta_e + 1500. / (theta_e * theta_e));
}

double fiducial_phi_hi(double theta_e, double theta_e_min, double k, double m) {
        return 37.14 * pow(sin((theta_e - theta_e_min) * 0.01745),
                           k + m / theta_e + 1500. / (theta_e * theta_e));
}

void MyProof::Begin(TTree * /*tree*/) {
        // The Begin() function is called at the start of the query.
        // When running with PROOF Begin() is only called on the client.
        // The tree argument is deprecated (on PROOF 0 is passed).

        TString option = GetOption();
        tBegin.Set();
        printf("*==* ---------- Begin of Job ----------");
        tBegin.Print();
}

void MyProof::SlaveBegin(TTree * /*tree*/) {
        // The SlaveBegin() function is called after the Begin() function.
        // When running with PROOF SlaveBegin() is called on each slave server.
        // The tree argument is deprecated (on PROOF 0 is passed).

        TString option = GetOption();
        printf("*==* ---------- Slave of Job ---------- \n");
        mom = new TH1D("mom", "mom", 500, 0, 5);
        fiducial = new TH2D("fiducial", "fiducial", 500, -PI, PI, 500, 0, 3);
        fiducial_cut =
                new TH2D("fiducial_cut", "fiducial_cut", 500, -PI, PI, 500, 0, 3);
        fOutput->AddAll(gDirectory->GetList());
}

Bool_t MyProof::Process(Long64_t entry) {
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

        // reset electron cut bool
        bool electron_cuts = true;
        // electron cuts
        electron_cuts &= (ec[0] > 0); // ``` ``` ``` ec
        electron_cuts &=
                ((int)id[0] == 11 || (int)id[0] == 0); // First particle is electron
        electron_cuts &=
                ((int)*gpart > 0); // Number of good particles is greater than 0
        electron_cuts &= ((int)stat[0] > 0); // First Particle hit stat
        electron_cuts &= ((int)q[0] == -1); // First particle is negative Q
        electron_cuts &= ((int)sc[0] > 0); // First Particle hit sc
        electron_cuts &= ((int)dc[0] > 0); // ``` ``` ``` dc
        electron_cuts &= ((int)dc_stat[dc[0] - 1] > 0);

        if (electron_cuts) {
                ++fNumberOfEvents;
                for (int i = 0; i < *gpart; i++) {
                        mom->Fill(p[i]);
                }

                e_mu_prime_3.SetXYZ(p[0] * cx[0], p[0] * cy[0], p[0] * cz[0]);
                e_mu_prime.SetVectM(e_mu_prime_3, 0.000511);
                double e_E = e_mu_prime.Energy();
                double costheta_e = e_mu_prime.CosTheta();
                double theta_e = e_mu_prime.Theta();
                double phi_e = e_mu_prime.Phi();
                double sector_width = PI / 3.0;
                int esector_index;
                for (esector_index = -3; (esector_index < 4); ++esector_index) {
                        if ((-0.5 * sector_width + esector_index * sector_width) <= phi_e &&
                            phi_e <= (0.5 * sector_width + esector_index * sector_width)) {
                                break;
                        }
                }

                int esector;
                if ((esector_index < 0)) {
                        esector = (esector_index + 7);
                } else {
                        esector = (esector_index + 1);
                }
                // phi angle relative to the sector:
                double phi_e_rel = phi_e - sector_width * esector_index;
                double efid_theta_min = 9.5 + 17.0 / (e_mu_prime.P() + 0.17);
                double efid_k = 0.705 + 1.1 * e_mu_prime.P();
                double efid_m = -63.5 + (-30.0 * e_mu_prime.P());
                double efid_phi_hi =
                        D2R * fiducial_phi_hi(theta_e / D2R, efid_theta_min, efid_k, efid_m);
                double efid_phi_lo =
                        D2R * fiducial_phi_lo(theta_e / D2R, efid_theta_min, efid_k, efid_m);

                bool efid_passes_cut = efid_phi_lo <= phi_e_rel && phi_e_rel <= efid_phi_hi;

                if (efid_passes_cut) {
                        fiducial_cut->Fill(phi_e, theta_e);
                } else {
                        fiducial->Fill(phi_e, theta_e);
                }
        }
        return kTRUE;
}

void MyProof::SlaveTerminate() {
        // The SlaveTerminate() function is called after all entries or objects
        // have been processed. When running with PROOF SlaveTerminate() is called
        // on each slave server.
        printf("\n *==* ---------- End of Slave Job ----------   ");
        tNow.Set();
        tNow.Print();
}

void MyProof::Terminate() {
        // The Terminate() function is the last function to be called during
        // a query. It always runs on the client, it can be used to present
        // the results graphically or save the results to file.

        TFile hfile("MyProof_Result.root", "RECREATE", "Results");
        fOutput->Write();

        tNow.Set();
        printf("*==* ---------- End of Job ---------- ");
        tNow.Print();
}
