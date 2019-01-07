//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jan  7 09:44:26 2019 by ROOT version 6.15/01
// from TTree h10/h10
// found on file: /Users/tylern/Data/e1d/new_cook/clas_022832.A01.root
//////////////////////////////////////////////////////////

#ifndef MySelector_h
#define MySelector_h

#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

// Headers needed by this particular selector

// Calcuating Q^2
//	Gotten from t channel
// -q^mu^2 = -(e^mu - e^mu')^2 = Q^2
double Q2_calc(TLorentzVector e_mu_prime) {
  TLorentzVector e_mu(0.0, 0.0, 4.81726, 4.81726);
  TLorentzVector q_mu = (e_mu - e_mu_prime);
  return -q_mu.Mag2();
}
//	Calcualting W
//	Gotten from s channel [(gamma + P)^2 == s == w^2]
//	Sqrt√[M_p^2 - Q^2 + 2 M_p gamma]
double W_calc(TLorentzVector e_mu_prime) {
  TLorentzVector e_mu(0.0, 0.0, 4.81726, 4.81726);
  TLorentzVector q_mu = (e_mu - e_mu_prime);
  TLorentzVector p_mu(0.0, 0.0, 0.0, 0.93827203);
  return (p_mu + q_mu).Mag();
}

class MySelector : public TSelector {
 public:
  TTreeReader fReader;  //! the tree reader
  TTree *fChain = 0;    //! pointer to the analyzed TTree or TChain
  TH2D *fWq2;
  TH1D *fW;

  // Readers to access the data (delete the ones you do not need).
  TTreeReaderValue<Int_t> npart = {fReader, "npart"};
  TTreeReaderValue<Int_t> evstat = {fReader, "evstat"};
  TTreeReaderValue<Int_t> intt = {fReader, "intt"};
  TTreeReaderValue<Int_t> evntid = {fReader, "evntid"};
  TTreeReaderValue<Int_t> evtype = {fReader, "evtype"};
  TTreeReaderValue<Int_t> evntclas = {fReader, "evntclas"};
  TTreeReaderValue<Int_t> evthel = {fReader, "evthel"};
  TTreeReaderValue<Int_t> evntclas2 = {fReader, "evntclas2"};
  TTreeReaderValue<Float_t> q_l = {fReader, "q_l"};
  TTreeReaderValue<Float_t> t_l = {fReader, "t_l"};
  TTreeReaderValue<Float_t> tr_time = {fReader, "tr_time"};
  TTreeReaderValue<Float_t> rf_time1 = {fReader, "rf_time1"};
  TTreeReaderValue<Float_t> rf_time2 = {fReader, "rf_time2"};
  TTreeReaderValue<Int_t> gpart = {fReader, "gpart"};
  TTreeReaderArray<Int_t> id = {fReader, "id"};
  TTreeReaderArray<Int_t> stat = {fReader, "stat"};
  TTreeReaderArray<Int_t> dc = {fReader, "dc"};
  TTreeReaderArray<Int_t> cc = {fReader, "cc"};
  TTreeReaderArray<Int_t> sc = {fReader, "sc"};
  TTreeReaderArray<Int_t> ec = {fReader, "ec"};
  TTreeReaderArray<Int_t> lec = {fReader, "lec"};
  TTreeReaderArray<Int_t> ccst = {fReader, "ccst"};
  TTreeReaderArray<Float_t> p = {fReader, "p"};
  TTreeReaderArray<Float_t> m = {fReader, "m"};
  TTreeReaderArray<Int_t> q = {fReader, "q"};
  TTreeReaderArray<Float_t> b = {fReader, "b"};
  TTreeReaderArray<Float_t> cx = {fReader, "cx"};
  TTreeReaderArray<Float_t> cy = {fReader, "cy"};
  TTreeReaderArray<Float_t> cz = {fReader, "cz"};
  TTreeReaderArray<Float_t> vx = {fReader, "vx"};
  TTreeReaderArray<Float_t> vy = {fReader, "vy"};
  TTreeReaderArray<Float_t> vz = {fReader, "vz"};
  TTreeReaderValue<Int_t> dc_part = {fReader, "dc_part"};
  TTreeReaderArray<Int_t> dc_sect = {fReader, "dc_sect"};
  TTreeReaderArray<Int_t> dc_trk = {fReader, "dc_trk"};
  TTreeReaderArray<Int_t> dc_stat = {fReader, "dc_stat"};
  TTreeReaderArray<Float_t> dc_vx = {fReader, "dc_vx"};
  TTreeReaderArray<Float_t> dc_vy = {fReader, "dc_vy"};
  TTreeReaderArray<Float_t> dc_vz = {fReader, "dc_vz"};
  TTreeReaderArray<Float_t> dc_vr = {fReader, "dc_vr"};
  TTreeReaderArray<Float_t> dc_xsc = {fReader, "dc_xsc"};
  TTreeReaderArray<Float_t> dc_ysc = {fReader, "dc_ysc"};
  TTreeReaderArray<Float_t> dc_zsc = {fReader, "dc_zsc"};
  TTreeReaderArray<Float_t> dc_cxsc = {fReader, "dc_cxsc"};
  TTreeReaderArray<Float_t> dc_cysc = {fReader, "dc_cysc"};
  TTreeReaderArray<Float_t> dc_czsc = {fReader, "dc_czsc"};
  TTreeReaderArray<Float_t> dc_c2 = {fReader, "dc_c2"};
  TTreeReaderValue<Int_t> ec_part = {fReader, "ec_part"};
  TTreeReaderArray<Int_t> ec_stat = {fReader, "ec_stat"};
  TTreeReaderArray<Int_t> ec_sect = {fReader, "ec_sect"};
  TTreeReaderArray<Int_t> ec_whol = {fReader, "ec_whol"};
  TTreeReaderArray<Int_t> ec_inst = {fReader, "ec_inst"};
  TTreeReaderArray<Int_t> ec_oust = {fReader, "ec_oust"};
  TTreeReaderArray<Float_t> etot = {fReader, "etot"};
  TTreeReaderArray<Float_t> ec_ei = {fReader, "ec_ei"};
  TTreeReaderArray<Float_t> ec_eo = {fReader, "ec_eo"};
  TTreeReaderArray<Float_t> ec_t = {fReader, "ec_t"};
  TTreeReaderArray<Float_t> ec_r = {fReader, "ec_r"};
  TTreeReaderArray<Float_t> ech_x = {fReader, "ech_x"};
  TTreeReaderArray<Float_t> ech_y = {fReader, "ech_y"};
  TTreeReaderArray<Float_t> ech_z = {fReader, "ech_z"};
  TTreeReaderArray<Float_t> ec_m2 = {fReader, "ec_m2"};
  TTreeReaderArray<Float_t> ec_m3 = {fReader, "ec_m3"};
  TTreeReaderArray<Float_t> ec_m4 = {fReader, "ec_m4"};
  TTreeReaderArray<Float_t> ec_c2 = {fReader, "ec_c2"};
  TTreeReaderValue<Int_t> sc_part = {fReader, "sc_part"};
  TTreeReaderArray<Int_t> sc_sect = {fReader, "sc_sect"};
  TTreeReaderArray<Int_t> sc_hit = {fReader, "sc_hit"};
  TTreeReaderArray<Int_t> sc_pd = {fReader, "sc_pd"};
  TTreeReaderArray<Int_t> sc_stat = {fReader, "sc_stat"};
  TTreeReaderArray<Float_t> edep = {fReader, "edep"};
  TTreeReaderArray<Float_t> sc_t = {fReader, "sc_t"};
  TTreeReaderArray<Float_t> sc_r = {fReader, "sc_r"};
  TTreeReaderArray<Float_t> sc_c2 = {fReader, "sc_c2"};
  TTreeReaderValue<Int_t> cc_part = {fReader, "cc_part"};
  TTreeReaderArray<Int_t> cc_sect = {fReader, "cc_sect"};
  TTreeReaderArray<Int_t> cc_hit = {fReader, "cc_hit"};
  TTreeReaderArray<Int_t> cc_segm = {fReader, "cc_segm"};
  TTreeReaderArray<Int_t> nphe = {fReader, "nphe"};
  TTreeReaderArray<Float_t> cc_t = {fReader, "cc_t"};
  TTreeReaderArray<Float_t> cc_r = {fReader, "cc_r"};
  TTreeReaderArray<Float_t> cc_c2 = {fReader, "cc_c2"};
  TTreeReaderValue<Int_t> lac_part = {fReader, "lac_part"};
  TTreeReaderArray<Int_t> lec_sect = {fReader, "lec_sect"};
  TTreeReaderArray<Int_t> lec_hit = {fReader, "lec_hit"};
  TTreeReaderArray<Int_t> lec_stat = {fReader, "lec_stat"};
  TTreeReaderArray<Float_t> lec_etot = {fReader, "lec_etot"};
  TTreeReaderArray<Float_t> lec_ein = {fReader, "lec_ein"};
  TTreeReaderArray<Float_t> lec_t = {fReader, "lec_t"};
  TTreeReaderArray<Float_t> lec_r = {fReader, "lec_r"};
  TTreeReaderArray<Float_t> lec_x = {fReader, "lec_x"};
  TTreeReaderArray<Float_t> lec_y = {fReader, "lec_y"};
  TTreeReaderArray<Float_t> lec_z = {fReader, "lec_z"};
  TTreeReaderArray<Float_t> lec_c2 = {fReader, "lec_c2"};
  TTreeReaderValue<Int_t> st_part = {fReader, "st_part"};
  TTreeReaderArray<Int_t> st_status = {fReader, "st_status"};
  TTreeReaderArray<Float_t> st_time = {fReader, "st_time"};
  TTreeReaderArray<Float_t> st_rtrk = {fReader, "st_rtrk"};

  MySelector(TTree * /*tree*/ = 0) {
    fWq2 = 0;
    fW = 0;
  }
  virtual ~MySelector() {}
  virtual Int_t Version() const { return 2; }
  virtual void Begin(TTree *tree);
  virtual void SlaveBegin(TTree *tree);
  virtual void Init(TTree *tree);
  virtual Bool_t Notify();
  virtual Bool_t Process(Long64_t entry);
  virtual Int_t GetEntry(Long64_t entry, Int_t getall = 0) {
    return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0;
  }
  virtual void SetOption(const char *option) { fOption = option; }
  virtual void SetObject(TObject *obj) { fObject = obj; }
  virtual void SetInputList(TList *input) { fInput = input; }
  virtual TList *GetOutputList() const { return fOutput; }
  virtual void SlaveTerminate();
  virtual void Terminate();

  ClassDef(MySelector, 0);
};

#endif

#ifdef MySelector_cxx
void MySelector::Init(TTree *tree) {
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the reader is initialized.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  fReader.SetTree(tree);
}

Bool_t MySelector::Notify() {
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

#endif  // #ifdef MySelector_cxx
