// A basic analysis program looking at
// a electron beam on a proton target
//
// Follow the TODO portions to get the analysis working properly
//
#include <TClonesArray.h>
#include <TObject.h>
#include "/usr/local/clas-software/analysis/ClasTool/include/TCCPBClass.h"
#include "/usr/local/clas-software/analysis/ClasTool/include/TDCPBClass.h"
#include "/usr/local/clas-software/analysis/ClasTool/include/TECPBClass.h"
#include "/usr/local/clas-software/analysis/ClasTool/include/THEADERClass.h"
#include "/usr/local/clas-software/analysis/ClasTool/include/TLCPBClass.h"
#include "/usr/local/clas-software/analysis/ClasTool/include/TSCPBClass.h"
#include "/usr/local/clas-software/analysis/ClasTool/include/TTGBIClass.h"
#include "/usr/local/clas-software/analysis/ClasTool/include/TVirtualData.h"

#define Square(x) ((x) * (x))

static const double BEAM = 4.81726;         // Beam energy in GeV
static const double MASS_P = 0.93827203;    // Mass in GeV
static const double MASS_E = 0.000511;      // Mass in GeV
static const double MASS_PIP = 0.13957018;  // Mass in GeV
const Int_t kMaxEVNT = 19;
const Int_t kMaxECPB = 18;
const Int_t kMaxSCPB = 18;
const Int_t kMaxDCPB = 18;
const Int_t kMaxCCPB = 7;
const Int_t kMaxLCPB = 11;
const Int_t kMaxTGBI = 1;

TLorentzVector e_mu(0.0, 0.0, TMath::Sqrt((BEAM * BEAM) - (MASS_E * MASS_E)), BEAM);
TLorentzVector p_mu(0.0, 0.0, 0.0, MASS_P);

double Q2_calc(TLorentzVector q_mu) { return -q_mu.Mag2(); }

double W_calc(TLorentzVector q_mu) { return (p_mu + q_mu).Mag(); }

// TCanvas c1("c1", "Plots", 1280, 720);

double Breit(double* x, double* par) { return par[2] * TMath::BreitWigner(x[0], par[0], par[1]); }

double missing_mass(TLorentzVector gamma_mu, TLorentzVector pip_mu) {
  TVector3 target_3;
  TLorentzVector target;
  // Set target vector
  target_3.SetXYZ(0.0, 0.0, 0.0);
  target.SetVectM(target_3, MASS_P);
  TLorentzVector reaction;
  reaction = (gamma_mu + target - pip_mu);

  return reaction.M();
}

void analyze() {
  const char* fout = "WvsQ2.root";
  TH1D* W_hist = new TH1D("W_hist", "W_hist", 500, 0.0, 3.0);
  TH1D* Q2_hist = new TH1D("Q2_hist", "Q2_hist", 500, 0.0, 4.0);
  TH2D* W_vs_Q2_hist = new TH2D("W_vs_Q2_hist", "W_vs_Q2_hist", 500, 0.0, 3.0, 500, 0.0, 4.0);
  TFile* OutputFile = new TFile(fout, "RECREATE");
  double W, Q2;
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(1111);
  gSystem->Load("/usr/local/clas-software/analysis/ClasTool/slib/Linux/libPartSieve.so");
  gSystem->Load("/usr/local/clas-software/analysis/ClasTool/slib/Linux/libFillBanks.so");
  gSystem->Load("/usr/local/clas-software/analysis/ClasTool/slib/Linux/libClasBanks.so");
  gSystem->Load("/usr/local/clas-software/analysis/ClasTool/slib/Linux/libDSTReader.so");
  gSystem->Load("/usr/local/clas-software/analysis/ClasTool/slib/Linux/libMapUtils.so");
  gSystem->Load("/usr/local/clas-software/analysis/ClasTool/slib/Linux/libFilter.so");
  gSystem->Load("/usr/local/clas-software/analysis/ClasTool/slib/Linux/libClasTool.so");
  gSystem->Load("/usr/local/clas-software/analysis/ClasTool/slib/Linux/libVirtualReader.so");

  TH1D* MM = new TH1D("mm", "mm", 500, 0, 3);

  TF1* bw = new TF1("bw", Breit, 0, 2, 3);
  /*
    Int_t NRun;
    Int_t NEvent;
    Int_t Time;
    Int_t Type;
    Int_t ROC;
    Int_t EvtClas;
    Int_t TrigBits;
    Int_t EStatus;
    Int_t TrgPrs;
    Int_t NPGP;
    Float_t FC;
    Float_t FCG;
    Float_t TG;
    Float_t STT;
    Float_t RF1;
    Int_t Latch1;
    Int_t Helicity_Scaler;
    Int_t Interrupt_Time;
    */
  Int_t EVNT_;
  UInt_t EVNT_fUniqueID[kMaxEVNT];  //[EVNT_]
  UInt_t EVNT_fBits[kMaxEVNT];      //[EVNT_]
  Int_t EVNT_Id[kMaxEVNT];          //[EVNT_]
  Char_t EVNT_Charge[kMaxEVNT];     //[EVNT_]
  Float_t EVNT_Betta[kMaxEVNT];     //[EVNT_]
  Float_t EVNT_Px[kMaxEVNT];        //[EVNT_]
  Float_t EVNT_Py[kMaxEVNT];        //[EVNT_]
  Float_t EVNT_Pz[kMaxEVNT];        //[EVNT_]
  Float_t EVNT_X[kMaxEVNT];         //[EVNT_]
  Float_t EVNT_Y[kMaxEVNT];         //[EVNT_]
  Float_t EVNT_Z[kMaxEVNT];         //[EVNT_]
  /*
  UChar_t EVNT_Dcstat[kMaxEVNT];    //[EVNT_]
  UChar_t EVNT_Ccstat[kMaxEVNT];    //[EVNT_]
  UChar_t EVNT_Scstat[kMaxEVNT];    //[EVNT_]
  UChar_t EVNT_Ecstat[kMaxEVNT];    //[EVNT_]
  UChar_t EVNT_Lcstat[kMaxEVNT];    //[EVNT_]
  UChar_t EVNT_Status[kMaxEVNT];    //[EVNT_]
  Int_t ECPB_;
  UInt_t ECPB_fUniqueID[kMaxECPB];  //[ECPB_]
  UInt_t ECPB_fBits[kMaxECPB];      //[ECPB_]
  Int_t ECPB_Scht[kMaxECPB];        //[ECPB_]
  Float_t ECPB_Etot[kMaxECPB];      //[ECPB_]
  Float_t ECPB_Ein[kMaxECPB];       //[ECPB_]
  Float_t ECPB_Eout[kMaxECPB];      //[ECPB_]
  Float_t ECPB_Time[kMaxECPB];      //[ECPB_]
  Float_t ECPB_Path[kMaxECPB];      //[ECPB_]
  Float_t ECPB_X[kMaxECPB];         //[ECPB_]
  Float_t ECPB_Y[kMaxECPB];         //[ECPB_]
  Float_t ECPB_Z[kMaxECPB];         //[ECPB_]
  Float_t ECPB_M2_hit[kMaxECPB];    //[ECPB_]
  Float_t ECPB_M3_hit[kMaxECPB];    //[ECPB_]
  Float_t ECPB_M4_hit[kMaxECPB];    //[ECPB_]
  Int_t ECPB_Innstr[kMaxECPB];      //[ECPB_]
  Int_t ECPB_Outstr[kMaxECPB];      //[ECPB_]
  Float_t ECPB_Chi2ec[kMaxECPB];    //[ECPB_]
  Int_t ECPB_Status[kMaxECPB];      //[ECPB_]
  Int_t SCPB_;
  UInt_t SCPB_fUniqueID[kMaxSCPB];  //[SCPB_]
  UInt_t SCPB_fBits[kMaxSCPB];      //[SCPB_]
  Int_t SCPB_Scpdht[kMaxSCPB];      //[SCPB_]
  Float_t SCPB_Edep[kMaxSCPB];      //[SCPB_]
  Float_t SCPB_Time[kMaxSCPB];      //[SCPB_]
  Float_t SCPB_Path[kMaxSCPB];      //[SCPB_]
  Float_t SCPB_Chi2sc[kMaxSCPB];    //[SCPB_]
  Int_t SCPB_Status[kMaxSCPB];      //[SCPB_]
  Int_t DCPB_;
  UInt_t DCPB_fUniqueID[kMaxDCPB];  //[DCPB_]
  UInt_t DCPB_fBits[kMaxDCPB];      //[DCPB_]
  Int_t DCPB_Sctr[kMaxDCPB];        //[DCPB_]
  Float_t DCPB_X_sc[kMaxDCPB];      //[DCPB_]
  Float_t DCPB_Y_sc[kMaxDCPB];      //[DCPB_]
  Float_t DCPB_Z_sc[kMaxDCPB];      //[DCPB_]
  Float_t DCPB_Cx_sc[kMaxDCPB];     //[DCPB_]
  Float_t DCPB_Cy_sc[kMaxDCPB];     //[DCPB_]
  Float_t DCPB_Cz_sc[kMaxDCPB];     //[DCPB_]
  Float_t DCPB_X_ec[kMaxDCPB];      //[DCPB_]
  Float_t DCPB_Y_ec[kMaxDCPB];      //[DCPB_]
  Float_t DCPB_Z_ec[kMaxDCPB];      //[DCPB_]
  Float_t DCPB_Th_cc[kMaxDCPB];     //[DCPB_]
  Float_t DCPB_Chi2[kMaxDCPB];      //[DCPB_]
  Int_t DCPB_Status[kMaxDCPB];      //[DCPB_]
  Int_t CCPB_;
  UInt_t CCPB_fUniqueID[kMaxCCPB];  //[CCPB_]
  UInt_t CCPB_fBits[kMaxCCPB];      //[CCPB_]
  Int_t CCPB_Scsght[kMaxCCPB];      //[CCPB_]
  Float_t CCPB_Nphe[kMaxCCPB];      //[CCPB_]
  Float_t CCPB_Time[kMaxCCPB];      //[CCPB_]
  Float_t CCPB_Path[kMaxCCPB];      //[CCPB_]
  Float_t CCPB_Chi2cc[kMaxCCPB];    //[CCPB_]
  Int_t CCPB_Status[kMaxCCPB];      //[CCPB_]
  Int_t LCPB_;
  UInt_t LCPB_fUniqueID[kMaxLCPB];  //[LCPB_]
  UInt_t LCPB_fBits[kMaxLCPB];      //[LCPB_]
  Int_t LCPB_Scht[kMaxLCPB];        //[LCPB_]
  Float_t LCPB_Etot[kMaxLCPB];      //[LCPB_]
  Float_t LCPB_Time[kMaxLCPB];      //[LCPB_]
  Float_t LCPB_Path[kMaxLCPB];      //[LCPB_]
  Float_t LCPB_X[kMaxLCPB];         //[LCPB_]
  Float_t LCPB_Y[kMaxLCPB];         //[LCPB_]
  Float_t LCPB_Z[kMaxLCPB];         //[LCPB_]
  Float_t LCPB_Chi2lc[kMaxLCPB];    //[LCPB_]
  Int_t LCPB_Status[kMaxLCPB];      //[LCPB_]
  Float_t LCPB_Ein[kMaxLCPB];       //[LCPB_]
  Int_t TGBI_;
  UInt_t TGBI_fUniqueID[kMaxTGBI];       //[TGBI_]
  UInt_t TGBI_fBits[kMaxTGBI];           //[TGBI_]
  Int_t TGBI_Latch1[kMaxTGBI];           //[TGBI_]
  Int_t TGBI_Helicity_scaler[kMaxTGBI];  //[TGBI_]
  Int_t TGBI_Interrupt_time[kMaxTGBI];   //[TGBI_]
  Int_t TGBI_Latch2[kMaxTGBI];           //[TGBI_]
  Int_t TGBI_Level3[kMaxTGBI];           //[TGBI_]
*/
  // Load chain from branch h10
  TChain event("CLASEVENT");
  event.Add("/mnt/500GB/rdst/rdst23040_07.root");
  /*
    event.SetBranchAddress("NRun", &NRun);
    event.SetBranchAddress("NEvent", &NEvent);
    event.SetBranchAddress("Time", &Time);
    event.SetBranchAddress("Type", &Type);
    event.SetBranchAddress("ROC", &ROC);
    event.SetBranchAddress("EvtClas", &EvtClas);
    event.SetBranchAddress("TrigBits", &TrigBits);
    event.SetBranchAddress("EStatus", &EStatus);
    event.SetBranchAddress("TrgPrs", &TrgPrs);
    event.SetBranchAddress("NPGP", &NPGP);
    event.SetBranchAddress("FC", &FC);
    event.SetBranchAddress("FCG", &FCG);
    event.SetBranchAddress("TG", &TG);
    event.SetBranchAddress("STT", &STT);
    event.SetBranchAddress("RF1", &RF1);
    event.SetBranchAddress("Latch1", &Latch1);
    event.SetBranchAddress("Helicity_Scaler", &Helicity_Scaler);
    event.SetBranchAddress("Interrupt_Time", &Interrupt_Time);
    */
  event.SetBranchAddress("EVNT", &EVNT_);
  event.SetBranchAddress("EVNT.fUniqueID", EVNT_fUniqueID);
  event.SetBranchAddress("EVNT.fBits", EVNT_fBits);
  event.SetBranchAddress("EVNT.Id", EVNT_Id);
  event.SetBranchAddress("EVNT.Charge", EVNT_Charge);
  event.SetBranchAddress("EVNT.Betta", EVNT_Betta);
  event.SetBranchAddress("EVNT.Px", EVNT_Px);
  event.SetBranchAddress("EVNT.Py", EVNT_Py);
  event.SetBranchAddress("EVNT.Pz", EVNT_Pz);
  event.SetBranchAddress("EVNT.X", EVNT_X);
  event.SetBranchAddress("EVNT.Y", EVNT_Y);
  event.SetBranchAddress("EVNT.Z", EVNT_Z);
  /*
    event.SetBranchAddress("EVNT.Dcstat", EVNT_Dcstat);
    event.SetBranchAddress("EVNT.Ccstat", EVNT_Ccstat);
    event.SetBranchAddress("EVNT.Scstat", EVNT_Scstat);
    event.SetBranchAddress("EVNT.Ecstat", EVNT_Ecstat);
    event.SetBranchAddress("EVNT.Lcstat", EVNT_Lcstat);
    event.SetBranchAddress("EVNT.Status", EVNT_Status);
    event.SetBranchAddress("ECPB", &ECPB_);
    event.SetBranchAddress("ECPB.fUniqueID", ECPB_fUniqueID);
    event.SetBranchAddress("ECPB.fBits", ECPB_fBits);
    event.SetBranchAddress("ECPB.Scht", ECPB_Scht);
    event.SetBranchAddress("ECPB.Etot", ECPB_Etot);
    event.SetBranchAddress("ECPB.Ein", ECPB_Ein);
    event.SetBranchAddress("ECPB.Eout", ECPB_Eout);
    event.SetBranchAddress("ECPB.Time", ECPB_Time);
    event.SetBranchAddress("ECPB.Path", ECPB_Path);
    event.SetBranchAddress("ECPB.X", ECPB_X);
    event.SetBranchAddress("ECPB.Y", ECPB_Y);
    event.SetBranchAddress("ECPB.Z", ECPB_Z);
    event.SetBranchAddress("ECPB.M2_hit", ECPB_M2_hit);
    event.SetBranchAddress("ECPB.M3_hit", ECPB_M3_hit);
    event.SetBranchAddress("ECPB.M4_hit", ECPB_M4_hit);
    event.SetBranchAddress("ECPB.Innstr", ECPB_Innstr);
    event.SetBranchAddress("ECPB.Outstr", ECPB_Outstr);
    event.SetBranchAddress("ECPB.Chi2ec", ECPB_Chi2ec);
    event.SetBranchAddress("ECPB.Status", ECPB_Status);
    event.SetBranchAddress("SCPB", &SCPB_);
    event.SetBranchAddress("SCPB.fUniqueID", SCPB_fUniqueID);
    event.SetBranchAddress("SCPB.fBits", SCPB_fBits);
    event.SetBranchAddress("SCPB.Scpdht", SCPB_Scpdht);
    event.SetBranchAddress("SCPB.Edep", SCPB_Edep);
    event.SetBranchAddress("SCPB.Time", SCPB_Time);
    event.SetBranchAddress("SCPB.Path", SCPB_Path);
    event.SetBranchAddress("SCPB.Chi2sc", SCPB_Chi2sc);
    event.SetBranchAddress("SCPB.Status", SCPB_Status);
    event.SetBranchAddress("DCPB", &DCPB_);
    event.SetBranchAddress("DCPB.fUniqueID", DCPB_fUniqueID);
    event.SetBranchAddress("DCPB.fBits", DCPB_fBits);
    event.SetBranchAddress("DCPB.Sctr", DCPB_Sctr);
    event.SetBranchAddress("DCPB.X_sc", DCPB_X_sc);
    event.SetBranchAddress("DCPB.Y_sc", DCPB_Y_sc);
    event.SetBranchAddress("DCPB.Z_sc", DCPB_Z_sc);
    event.SetBranchAddress("DCPB.Cx_sc", DCPB_Cx_sc);
    event.SetBranchAddress("DCPB.Cy_sc", DCPB_Cy_sc);
    event.SetBranchAddress("DCPB.Cz_sc", DCPB_Cz_sc);
    event.SetBranchAddress("DCPB.X_ec", DCPB_X_ec);
    event.SetBranchAddress("DCPB.Y_ec", DCPB_Y_ec);
    event.SetBranchAddress("DCPB.Z_ec", DCPB_Z_ec);
    event.SetBranchAddress("DCPB.Th_cc", DCPB_Th_cc);
    event.SetBranchAddress("DCPB.Chi2", DCPB_Chi2);
    event.SetBranchAddress("DCPB.Status", DCPB_Status);
    event.SetBranchAddress("CCPB", &CCPB_);
    event.SetBranchAddress("CCPB.fUniqueID", CCPB_fUniqueID);
    event.SetBranchAddress("CCPB.fBits", CCPB_fBits);
    event.SetBranchAddress("CCPB.Scsght", CCPB_Scsght);
    event.SetBranchAddress("CCPB.Nphe", CCPB_Nphe);
    event.SetBranchAddress("CCPB.Time", CCPB_Time);
    event.SetBranchAddress("CCPB.Path", CCPB_Path);
    event.SetBranchAddress("CCPB.Chi2cc", CCPB_Chi2cc);
    event.SetBranchAddress("CCPB.Status", CCPB_Status);
    event.SetBranchAddress("LCPB", &LCPB_);
    event.SetBranchAddress("LCPB.fUniqueID", LCPB_fUniqueID);
    event.SetBranchAddress("LCPB.fBits", LCPB_fBits);
    event.SetBranchAddress("LCPB.Scht", LCPB_Scht);
    event.SetBranchAddress("LCPB.Etot", LCPB_Etot);
    event.SetBranchAddress("LCPB.Time", LCPB_Time);
    event.SetBranchAddress("LCPB.Path", LCPB_Path);
    event.SetBranchAddress("LCPB.X", LCPB_X);
    event.SetBranchAddress("LCPB.Y", LCPB_Y);
    event.SetBranchAddress("LCPB.Z", LCPB_Z);
    event.SetBranchAddress("LCPB.Chi2lc", LCPB_Chi2lc);
    event.SetBranchAddress("LCPB.Status", LCPB_Status);
    event.SetBranchAddress("LCPB.Ein", LCPB_Ein);
    event.SetBranchAddress("TGBI", &TGBI_);
    event.SetBranchAddress("TGBI.fUniqueID", TGBI_fUniqueID);
    event.SetBranchAddress("TGBI.fBits", TGBI_fBits);
    event.SetBranchAddress("TGBI.Latch1", TGBI_Latch1);
    event.SetBranchAddress("TGBI.Helicity_scaler", TGBI_Helicity_scaler);
    event.SetBranchAddress("TGBI.Interrupt_time", TGBI_Interrupt_time);
    event.SetBranchAddress("TGBI.Latch2", TGBI_Latch2);
    event.SetBranchAddress("TGBI.Level3", TGBI_Level3);
  */
  int num_of_events = (int)event.GetEntries();
  cout << num_of_events << endl;
  for (int current_event = 0; current_event < num_of_events; current_event++) {
    if (current_event > 10) break;
    event.GetEntry(current_event);

    // Setup scattered electron 4 vector
    TVector3 e_mu_prime_3;
    TLorentzVector e_mu_prime;
    TLorentzVector e_mu(0.0, 0.0, sqrt(Square(BEAM) - Square(MASS_E)), BEAM);
    for (int i = 0; i < EVNT_; i++) {
      cout << EVNT_Id[i] << endl;
    }

    if (EVNT_Id[0] == 0) {
      e_mu_prime_3.SetXYZ(EVNT_Px[0], EVNT_Py[0], EVNT_Pz[0]);
      e_mu_prime.SetVectM(e_mu_prime_3, MASS_E);

      q_mu = (e_mu - e_mu_prime);
      Q2 = Q2_calc(q_mu);
      W = W_calc(q_mu);
      cout << Q2 << " , " << W << endl;
      W_hist->Fill(W);
      Q2_hist->Fill(Q2);
      W_vs_Q2_hist->Fill(W, Q2);
    }
  }
  //
  // end stuff
  event.Reset();
  OutputFile->cd();
  W_hist->Write();
  Q2_hist->Write();
  W_vs_Q2_hist->Write();
  OutputFile->Close();
}
