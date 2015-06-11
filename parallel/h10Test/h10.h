//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jun 10 15:48:16 2011 by ROOT version 5.31/01
// from TChain h10/
//////////////////////////////////////////////////////////

#ifndef h10_h
#define h10_h

#include <iostream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH1.h>
#include <TH2.h>

using namespace std;

class h10 : public TSelector {
public :

   /* **********************************************************************
      ************************* FIRST PASS SUMMARY ITEMS *******************
      ********************************************************************** */
   Int_t filenumber;
   Int_t nfc;
   Double_t ql[200];
   Double_t ql_sum;
   Int_t nSurviveE, nSurviveEP, nSurviveEPPip, nSurviveEPPim, nSurviveEPPipPim;
   Int_t nSurviveE_per, nSurviveEP_per, nSurviveEPPip_per, nSurviveEPPim_per, nSurviveEPPipPim_per;
   class TH1I *h_nParts, *h_nPos, *h_nNeg, *h_nE, *h_nP, *h_nPip, *h_nPim, *h_nSurvive;
   class TH1F *h_nSurviveE_per, *h_nSurviveEP_per, *h_nSurviveEPPip_per, *h_nSurviveEPPim_per, *h_nSurviveEPPipPim_per;
   TH2F *h_bvp[4];
   /* ********************************************************************** */

   char _esf_p[4][6][128], _esf_inout[4][6][128], _nph[4][6][128];

   Float_t _xlow_h_eSF_p, _xhigh_h_eSF_p;
   Int_t _xbins_h_eSF_p;
   Float_t _ylow_h_eSF_p, _yhigh_h_eSF_p;
   Int_t _ybins_h_eSF_p;
   char _name_h_eSF_p[100];
   char _title_h_eSF_p[100];
   char _xAxisTitle_h_eSF_p[100];
   char _yAxisTitle_h_eSF_p[100];

   Float_t _xlow_h_eSFout_eSFin, _xhigh_h_eSFout_eSFin;
   Int_t _xbins_h_eSFout_eSFin;
   Float_t _ylow_h_eSFout_eSFin, _yhigh_h_eSFout_eSFin;
   Int_t _ybins_h_eSFout_eSFin;
   char _name_h_eSFout_eSFin[100];
   char _title_h_eSFout_eSFin[100];
   char _xAxisTitle_h_eSFout_eSFin[100];
   char _yAxisTitle_h_eSFout_eSFin[100];

   Float_t _xlow_h_nPhe, _xhigh_h_nPhe;
   Int_t _xbins_h_nPhe;
   char _name_h_nPhe[100];
   char _title_h_nPhe[100];
   char _xAxisTitle_h_nPhe[100];
   char _yAxisTitle_h_nPhe[100];
   TH2F* _h_eSF_p[4][6];
   TH2F* _h_eSFout_eSFin[4][6];
   TH1F* _h_nPhe[4][6];

   void InitEid();
   /* ********************************* */

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   UChar_t         npart;
   UChar_t         evstat;
   UInt_t          evntid;
   Char_t          evtype;
   Char_t          evntclas;
   Char_t          evthel;
   Int_t           evntclas2;
   Float_t         q_l;
   Float_t         t_l;
   Float_t         tr_time;
   Float_t         rf_time1;
   Float_t         rf_time2;
   Int_t           gpart;
   Short_t         id[26];   //[gpart]
   Char_t          stat[26];   //[gpart]
   UChar_t         dc[26];   //[gpart]
   UChar_t         cc[26];   //[gpart]
   UChar_t         sc[26];   //[gpart]
   UChar_t         ec[26];   //[gpart]
   UChar_t         lec[26];   //[gpart]
   Float_t         p[26];   //[gpart]
   Float_t         m[26];   //[gpart]
   Char_t          q[26];   //[gpart]
   Float_t         beta[26];   //[gpart]
   Float_t         cx[26];   //[gpart]
   Float_t         cy[26];   //[gpart]
   Float_t         cz[26];   //[gpart]
   Float_t         vx[26];   //[gpart]
   Float_t         vy[26];   //[gpart]
   Float_t         vz[26];   //[gpart]
   Int_t           dc_part;
   UChar_t         dc_sect[20];   //[dc_part]
   UChar_t         dc_trk[20];   //[dc_part]
   Char_t          dc_stat[20];   //[dc_part]
   Float_t         dc_xsc[20];   //[dc_part]
   Float_t         dc_ysc[20];   //[dc_part]
   Float_t         dc_zsc[20];   //[dc_part]
   Float_t         dc_cxsc[20];   //[dc_part]
   Float_t         dc_cysc[20];   //[dc_part]
   Float_t         dc_czsc[20];   //[dc_part]
   Float_t         dc_xec[20];   //[dc_part]
   Float_t         dc_yec[20];   //[dc_part]
   Float_t         dc_zec[20];   //[dc_part]
   Float_t         dc_thcc[20];   //[dc_part]
   Float_t         dc_c2[20];   //[dc_part]
   Int_t           ec_part;
   UShort_t        ec_stat[22];   //[ec_part]
   UChar_t         ec_sect[22];   //[ec_part]
   Int_t           ec_whol[22];   //[ec_part]
   Int_t           ec_inst[22];   //[ec_part]
   Int_t           ec_oust[22];   //[ec_part]
   Float_t         etot[22];   //[ec_part]
   Float_t         ec_ei[22];   //[ec_part]
   Float_t         ec_eo[22];   //[ec_part]
   Float_t         ec_t[22];   //[ec_part]
   Float_t         ec_r[22];   //[ec_part]
   Float_t         ech_x[22];   //[ec_part]
   Float_t         ech_y[22];   //[ec_part]
   Float_t         ech_z[22];   //[ec_part]
   Float_t         ec_m2[22];   //[ec_part]
   Float_t         ec_m3[22];   //[ec_part]
   Float_t         ec_m4[22];   //[ec_part]
   Float_t         ec_c2[22];   //[ec_part]
   Int_t           sc_part;
   UChar_t         sc_sect[22];   //[sc_part]
   UChar_t         sc_hit[22];   //[sc_part]
   UChar_t         sc_pd[22];   //[sc_part]
   UChar_t         sc_stat[22];   //[sc_part]
   Float_t         edep[22];   //[sc_part]
   Float_t         sc_t[22];   //[sc_part]
   Float_t         sc_r[22];   //[sc_part]
   Float_t         sc_c2[22];   //[sc_part]
   Int_t           cc_part;
   UChar_t         cc_sect[20];   //[cc_part]
   UChar_t         cc_hit[20];   //[cc_part]
   Int_t           cc_segm[20];   //[cc_part]
   UShort_t        nphe[20];   //[cc_part]
   Float_t         cc_t[20];   //[cc_part]
   Float_t         cc_r[20];   //[cc_part]
   Float_t         cc_c2[20];   //[cc_part]
   Int_t           lac_part;
   Int_t           lec_sect[20];   //[lac_part]
   Int_t           lec_hit[20];   //[lac_part]
   Int_t           lec_stat[20];   //[lac_part]
   Float_t         lec_etot[20];   //[lac_part]
   Float_t         lec_t[20];   //[lac_part]
   Float_t         lec_r[20];   //[lac_part]
   Float_t         lec_x[20];   //[lac_part]
   Float_t         lec_y[20];   //[lac_part]
   Float_t         lec_z[20];   //[lac_part]
   Float_t         lec_c2[20];   //[lac_part]

   // List of branches
   TBranch        *b_npart;   //!
   TBranch        *b_evstat;   //!
   TBranch        *b_evntid;   //!
   TBranch        *b_evtype;   //!
   TBranch        *b_evntclas;   //!
   TBranch        *b_evthel;   //!
   TBranch        *b_evntclas2;   //!
   TBranch        *b_q_l;   //!
   TBranch        *b_t_l;   //!
   TBranch        *b_tr_time;   //!
   TBranch        *b_rf_time1;   //!
   TBranch        *b_rf_time2;   //!
   TBranch        *b_gpart;   //!
   TBranch        *b_id;   //!
   TBranch        *b_stat;   //!
   TBranch        *b_dc;   //!
   TBranch        *b_cc;   //!
   TBranch        *b_sc;   //!
   TBranch        *b_ec;   //!
   TBranch        *b_lec;   //!
   TBranch        *b_p;   //!
   TBranch        *b_m;   //!
   TBranch        *b_q;   //!
   TBranch        *b_beta;   //!
   TBranch        *b_cx;   //!
   TBranch        *b_cy;   //!
   TBranch        *b_cz;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_dc_part;   //!
   TBranch        *b_dc_sect;   //!
   TBranch        *b_dc_trk;   //!
   TBranch        *b_dc_stat;   //!
   TBranch        *b_dc_xsc;   //!
   TBranch        *b_dc_ysc;   //!
   TBranch        *b_dc_zsc;   //!
   TBranch        *b_dc_cxsc;   //!
   TBranch        *b_dc_cysc;   //!
   TBranch        *b_dc_czsc;   //!
   TBranch        *b_dc_xec;   //!
   TBranch        *b_dc_yec;   //!
   TBranch        *b_dc_zec;   //!
   TBranch        *b_dc_thcc;   //!
   TBranch        *b_dc_c2;   //!
   TBranch        *b_ec_part;   //!
   TBranch        *b_ec_stat;   //!
   TBranch        *b_ec_sect;   //!
   TBranch        *b_ec_whol;   //!
   TBranch        *b_ec_inst;   //!
   TBranch        *b_ec_oust;   //!
   TBranch        *b_etot;   //!
   TBranch        *b_ec_ei;   //!
   TBranch        *b_ec_eo;   //!
   TBranch        *b_ec_t;   //!
   TBranch        *b_ec_r;   //!
   TBranch        *b_ech_x;   //!
   TBranch        *b_ech_y;   //!
   TBranch        *b_ech_z;   //!
   TBranch        *b_ec_m2;   //!
   TBranch        *b_ec_m3;   //!
   TBranch        *b_ec_m4;   //!
   TBranch        *b_ec_c2;   //!
   TBranch        *b_sc_part;   //!
   TBranch        *b_sc_sect;   //!
   TBranch        *b_sc_hit;   //!
   TBranch        *b_sc_pd;   //!
   TBranch        *b_sc_stat;   //!
   TBranch        *b_edep;   //!
   TBranch        *b_sc_t;   //!
   TBranch        *b_sc_r;   //!
   TBranch        *b_sc_c2;   //!
   TBranch        *b_cc_part;   //!
   TBranch        *b_cc_sect;   //!
   TBranch        *b_cc_hit;   //!
   TBranch        *b_cc_segm;   //!
   TBranch        *b_nphe;   //!
   TBranch        *b_cc_t;   //!
   TBranch        *b_cc_r;   //!
   TBranch        *b_cc_c2;   //!
   TBranch        *b_lac_part;   //!
   TBranch        *b_lec_sect;   //!
   TBranch        *b_lec_hit;   //!
   TBranch        *b_lec_stat;   //!
   TBranch        *b_lec_etot;   //!
   TBranch        *b_lec_t;   //!
   TBranch        *b_lec_r;   //!
   TBranch        *b_lec_x;   //!
   TBranch        *b_lec_y;   //!
   TBranch        *b_lec_z;   //!
   TBranch        *b_lec_c2;   //!

   h10(TTree * /*tree*/ =0) { }
   virtual ~h10() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(h10,0);
};

#endif

#ifdef h10_cxx
void h10::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).
   cout << "Init" << endl;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   //if (fChain->GetBranch("evthel")) {
   if (fChain->GetBranch("h10")) {
      fChain->SetBranchAddress("evstat", &evstat, &b_evstat);
      fChain->SetBranchAddress("evthel", &evthel, &b_evthel);
      fChain->SetBranchAddress("evntclas2", &evntclas2, &b_evntclas2);
      fChain->SetBranchAddress("rf_time2", &rf_time2, &b_rf_time2);
      fChain->SetBranchAddress("m", m, &b_m);
      fChain->SetBranchAddress("dc_xec", dc_xec, &b_dc_xec);
      fChain->SetBranchAddress("dc_yec", dc_yec, &b_dc_yec);
      fChain->SetBranchAddress("dc_zec", dc_zec, &b_dc_zec);
      fChain->SetBranchAddress("dc_thcc", dc_thcc, &b_dc_thcc);
      fChain->SetBranchAddress("ec_inst", ec_inst, &b_ec_inst);
      fChain->SetBranchAddress("ec_oust", ec_oust, &b_ec_oust);
      fChain->SetBranchAddress("ec_m2", ec_m2, &b_ec_m2);
      fChain->SetBranchAddress("ec_m3", ec_m3, &b_ec_m3);
      fChain->SetBranchAddress("ec_m4", ec_m4, &b_ec_m4);
      fChain->SetBranchAddress("lac_part", &lac_part, &b_lac_part);
      fChain->SetBranchAddress("lec_sect", lec_sect, &b_lec_sect);
      fChain->SetBranchAddress("lec_hit", lec_hit, &b_lec_hit);
      fChain->SetBranchAddress("lec_stat", lec_stat, &b_lec_stat);
      fChain->SetBranchAddress("lec_etot", lec_etot, &b_lec_etot);
      fChain->SetBranchAddress("lec_t", lec_t, &b_lec_t);
      fChain->SetBranchAddress("lec_r", lec_r, &b_lec_r);
      fChain->SetBranchAddress("lec_x", lec_x, &b_lec_x);
      fChain->SetBranchAddress("lec_y", lec_y, &b_lec_y);
      fChain->SetBranchAddress("lec_z", lec_z, &b_lec_z);
      fChain->SetBranchAddress("lec_c2", lec_c2, &b_lec_c2);
   }
   fChain->SetBranchAddress("npart", &npart, &b_npart);
   fChain->SetBranchAddress("evntid", &evntid, &b_evntid);
   fChain->SetBranchAddress("evtype", &evtype, &b_evtype);
   fChain->SetBranchAddress("evntclas", &evntclas, &b_evntclas);
   fChain->SetBranchAddress("q_l", &q_l, &b_q_l);
   fChain->SetBranchAddress("t_l", &t_l, &b_t_l);
   fChain->SetBranchAddress("tr_time", &tr_time, &b_tr_time);
   fChain->SetBranchAddress("rf_time1", &rf_time1, &b_rf_time1);
   fChain->SetBranchAddress("gpart", &gpart, &b_gpart);
   fChain->SetBranchAddress("id", id, &b_id);
   fChain->SetBranchAddress("stat", stat, &b_stat);
   fChain->SetBranchAddress("dc", dc, &b_dc);
   fChain->SetBranchAddress("cc", cc, &b_cc);
   fChain->SetBranchAddress("sc", sc, &b_sc);
   fChain->SetBranchAddress("ec", ec, &b_ec);
   fChain->SetBranchAddress("lec", lec, &b_lec);
   fChain->SetBranchAddress("p", p, &b_p);
   fChain->SetBranchAddress("q", q, &b_q);
   fChain->SetBranchAddress("b", beta, &b_beta);
   fChain->SetBranchAddress("cx", cx, &b_cx);
   fChain->SetBranchAddress("cy", cy, &b_cy);
   fChain->SetBranchAddress("cz", cz, &b_cz);
   fChain->SetBranchAddress("vx", vx, &b_vx);
   fChain->SetBranchAddress("vy", vy, &b_vy);
   fChain->SetBranchAddress("vz", vz, &b_vz);
   fChain->SetBranchAddress("dc_part", &dc_part, &b_dc_part);
   fChain->SetBranchAddress("dc_sect", dc_sect, &b_dc_sect);
   fChain->SetBranchAddress("dc_trk", dc_trk, &b_dc_trk);
   fChain->SetBranchAddress("dc_stat", dc_stat, &b_dc_stat);
   fChain->SetBranchAddress("dc_xsc", dc_xsc, &b_dc_xsc);
   fChain->SetBranchAddress("dc_ysc", dc_ysc, &b_dc_ysc);
   fChain->SetBranchAddress("dc_zsc", dc_zsc, &b_dc_zsc);
   fChain->SetBranchAddress("dc_cxsc", dc_cxsc, &b_dc_cxsc);
   fChain->SetBranchAddress("dc_cysc", dc_cysc, &b_dc_cysc);
   fChain->SetBranchAddress("dc_czsc", dc_czsc, &b_dc_czsc);
   fChain->SetBranchAddress("dc_c2", dc_c2, &b_dc_c2);
   fChain->SetBranchAddress("ec_part", &ec_part, &b_ec_part);
   fChain->SetBranchAddress("ec_stat", ec_stat, &b_ec_stat);
   fChain->SetBranchAddress("ec_sect", ec_sect, &b_ec_sect);
   fChain->SetBranchAddress("ec_whol", ec_whol, &b_ec_whol);
   fChain->SetBranchAddress("etot", etot, &b_etot);
   fChain->SetBranchAddress("ec_ei", ec_ei, &b_ec_ei);
   fChain->SetBranchAddress("ec_eo", ec_eo, &b_ec_eo);
   fChain->SetBranchAddress("ec_t", ec_t, &b_ec_t);
   fChain->SetBranchAddress("ec_r", ec_r, &b_ec_r);
   fChain->SetBranchAddress("ech_x", ech_x, &b_ech_x);
   fChain->SetBranchAddress("ech_y", ech_y, &b_ech_y);
   fChain->SetBranchAddress("ech_z", ech_z, &b_ech_z);
   fChain->SetBranchAddress("ec_c2", ec_c2, &b_ec_c2);
   fChain->SetBranchAddress("sc_part", &sc_part, &b_sc_part);
   fChain->SetBranchAddress("sc_sect", sc_sect, &b_sc_sect);
   fChain->SetBranchAddress("sc_hit", sc_hit, &b_sc_hit);
   fChain->SetBranchAddress("sc_pd", sc_pd, &b_sc_pd);
   fChain->SetBranchAddress("sc_stat", sc_stat, &b_sc_stat);
   fChain->SetBranchAddress("edep", edep, &b_edep);
   fChain->SetBranchAddress("sc_t", sc_t, &b_sc_t);
   fChain->SetBranchAddress("sc_r", sc_r, &b_sc_r);
   fChain->SetBranchAddress("sc_c2", sc_c2, &b_sc_c2);
   fChain->SetBranchAddress("cc_part", &cc_part, &b_cc_part);
   fChain->SetBranchAddress("cc_sect", cc_sect, &b_cc_sect);
   fChain->SetBranchAddress("cc_hit", cc_hit, &b_cc_hit);
   fChain->SetBranchAddress("cc_segm", cc_segm, &b_cc_segm);
   fChain->SetBranchAddress("nphe", nphe, &b_nphe);
   fChain->SetBranchAddress("cc_t", cc_t, &b_cc_t);
   fChain->SetBranchAddress("cc_r", cc_r, &b_cc_r);
   fChain->SetBranchAddress("cc_c2", cc_c2, &b_cc_c2);
}

Bool_t h10::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.
   cout << "Notify" << ++filenumber << endl;
   nfc = 0;
   ql_sum = 0;
   nSurviveE_per = nSurviveEP_per = nSurviveEPPip_per = nSurviveEPPim_per = nSurviveEPPipPim_per = 0;
   return kTRUE;
}

#endif // #ifdef h10_cxx