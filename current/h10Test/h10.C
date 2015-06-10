#define h10_cxx
// The class definition in h10.h has been generated automatically
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
// Root > T->Process("h10.C")
// Root > T->Process("h10.C","some options")
// Root > T->Process("h10.C+")
//

#include "h10.h"
#include <TBrowser.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TCanvas.h>

static const char all[] = "all";
static const char n1e[] = "e";
static const char n2pos1neg[] = "2+1-";
static const char n1pos2neg[] = "1+2-";
static const char n2pos2neg[] = "2+2-";
static const char n1e1p[] = "ep";
static const char n1e1p1pip[] = "eppi+";
static const char n1e1p1pim[] = "eppi-";
static const char n1e1p1pip1pim[] = "eppi+pi-";

void h10::InitEid() {

   _xlow_h_eSF_p = 0;
   _xhigh_h_eSF_p = 6;
   _ylow_h_eSF_p = 0;
   _yhigh_h_eSF_p = 0.7;
   _xbins_h_eSF_p = 100;
   _ybins_h_eSF_p = 50;
   sprintf(_name_h_eSF_p, "eSF_p");
   sprintf(_title_h_eSF_p, "e^{-} EC Sampling Fraction vs Momentum");
   sprintf(_xAxisTitle_h_eSF_p,"P (GeV)");
   sprintf(_yAxisTitle_h_eSF_p, "E_{Tot}/P");

   _xlow_h_eSFout_eSFin = 0;
   _xhigh_h_eSFout_eSFin = 0.4;
   _ylow_h_eSFout_eSFin = 0;
   _yhigh_h_eSFout_eSFin = 0.25;
   _xbins_h_eSFout_eSFin = 100;
   _ybins_h_eSFout_eSFin = 100;
   sprintf(_name_h_eSFout_eSFin, "eSFout_eSFin");
   sprintf(_title_h_eSFout_eSFin, "e^{-} EC Sampling Fraction: Outer vs Inner");
   sprintf(_xAxisTitle_h_eSFout_eSFin, "E_{in}/P");
   sprintf(_yAxisTitle_h_eSFout_eSFin, "E_{in}/P");

   _xlow_h_nPhe = 0;
   _xhigh_h_nPhe = 300;
   _xbins_h_nPhe = 300;
   sprintf(_name_h_nPhe, "nPhe");
   sprintf(_title_h_nPhe, "nPhe");
   sprintf(_xAxisTitle_h_nPhe, "nPhe*10");

   for(int i = 0; i < 6; i++){
      sprintf(_esf_p[0][i],"%s_s_%d",_name_h_eSF_p,i+1);
      sprintf(_esf_inout[0][i],"%s_s_%d",_name_h_eSFout_eSFin,i+1);
      sprintf(_nph[0][i],"%s_s_%d",_name_h_nPhe,i+1);
      sprintf(_esf_p[1][i],"%s_s_%d_%s",_name_h_eSF_p,i+1,"eOnly");
      sprintf(_esf_inout[1][i],"%s_s_%d_%s",_name_h_eSFout_eSFin,i+1,"eOnly");
      sprintf(_nph[1][i],"%s_s_%d_%s",_name_h_nPhe,i+1,"eOnly");
      sprintf(_esf_p[2][i],"%s_s_%d_%s",_name_h_eSF_p,i+1,"piOnly");
      sprintf(_esf_inout[2][i],"%s_s_%d_%s",_name_h_eSFout_eSFin,i+1,"piOnly");
      sprintf(_nph[2][i],"%s_s_%d_%s",_name_h_nPhe,i+1,"piOnly");
      sprintf(_esf_p[3][i],"%s_s_%d_%s",_name_h_eSF_p,i+1,"eOnly_not");
      sprintf(_esf_inout[3][i],"%s_s_%d_%s",_name_h_eSFout_eSFin,i+1,"eOnly_not");
      sprintf(_nph[3][i],"%s_s_%d_%s",_name_h_nPhe,i+1,"eOnly_not");
      // cout << _esf_p[0][i] << endl;
      // cout << _esf_inout[0][i] << endl;
      // cout << _nph[0][i] << endl;
   }
}

void h10::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   cout << "Begin" << endl;
   InitEid();
}

void h10::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   cout << "SlaveBegin" << endl;

   nSurviveE = nSurviveEP = nSurviveEPPip = nSurviveEPPim = nSurviveEPPipPim = 0;
   h_nSurviveE_per = new TH1F("h_nEperQ","e^{-} norm q",2000,0,100000);

   h_nParts = new TH1I("h_nParts", "# Events vs # Parts", 41,-0.5,40.5);
   h_nPos = new TH1I("h_nPos", "# Events vs # Positive Parts", 41,-0.5,40.5);
   h_nNeg = new TH1I("h_nNeg", "# Events vs # Negative Parts", 41,-0.5,40.5);
   h_nE = new TH1I("h_nE", "# Events vs # e^{-}", 41,-0.5,40.5);
   h_nP = new TH1I("h_nP", "# Events vs # p", 41,-0.5,40.5);
   h_nPip = new TH1I("h_nPip", "# Events vs # #pi^{+}", 41,-0.5,40.5);
   h_nPim = new TH1I("h_nPim", "# Events vs # #pi^{-}", 41,-0.5,40.5);
   h_nSurvive = new TH1I("h_nSurvive", "Particles per Event", 10,0.5,10.5);
   TAxis *xax = h_nSurvive->GetXaxis();
   xax->SetBinLabel(1,all);
   xax->SetBinLabel(2,n1e);
   xax->SetBinLabel(3,n2pos1neg);
   xax->SetBinLabel(4,n1pos2neg);
   xax->SetBinLabel(5,n1e1p);
   xax->SetBinLabel(6,n1e1p1pip);
   xax->SetBinLabel(7,n1e1p1pim);
   xax->SetBinLabel(8,n2pos2neg);
   xax->SetBinLabel(9,n1e1p1pip1pim);

   h_bvp[0] = new TH2F("hbvp","hbvp",510,-0.1,5,240,0,1.2);
   h_bvp[1] = new TH2F("hbvp_neg","hbvp negative",510,-0.1,5,240,0,1.2);
   h_bvp[2] = new TH2F("hbvp_pim","hbvp #pi^{-}",510,-0.1,5,240,0,1.2);
   h_bvp[3] = new TH2F("hbvp_pimNot","hbvp not #pi^{-}",510,-0.1,5,240,0,1.2);
   for (int i = 0; i < 4; i++) fOutput->Add(h_bvp[i]);
   /* ****************************************************** */
   InitEid();
   char esf_p_t[256];
   char esf_inout_t[256];
   char nph_t[256];
   for(int i = 0; i < 6; i++){
      //      sprintf(esf_p,"%s_s_%d",_name_h_eSF_p,i+1);
      sprintf(esf_p_t,"%s_s_%d",_title_h_eSF_p,i+1);
      //      sprintf(esf_inout,"%s_s_%d",_name_h_eSFout_eSFin,i+1);
      sprintf(esf_inout_t,"%s_s_%d",_title_h_eSFout_eSFin,i+1);
      //      sprintf(nph,"%s_s_%d",_name_h_nPhe,i+1);
      sprintf(nph_t,"%s_s_%d",_title_h_nPhe,i+1);
      //cout << "In sector" << i + 1 << endl;
      //cout << _esf_p[i] << endl;
      for (int j = 0; j < 4; j++) {
	 _h_eSF_p[j][i] = new TH2F(_esf_p[j][i], esf_p_t, _xbins_h_eSF_p, _xlow_h_eSF_p, _xhigh_h_eSF_p, _ybins_h_eSF_p, _ylow_h_eSF_p, _yhigh_h_eSF_p);
	 //cout << _esf_inout[i] << endl;
	 _h_eSFout_eSFin[j][i] = new TH2F(_esf_inout[j][i], esf_inout_t, _xbins_h_eSFout_eSFin, _xlow_h_eSFout_eSFin, _xhigh_h_eSFout_eSFin, _xbins_h_eSFout_eSFin, _xlow_h_eSFout_eSFin, _xhigh_h_eSFout_eSFin);
	 //cout << _nph[i] << endl;
	 _h_nPhe[j][i] = new TH1F(_nph[j][i], nph_t, _xbins_h_nPhe, _xlow_h_nPhe, _xhigh_h_nPhe);
	 fOutput->Add(_h_eSF_p[j][i]);
	 fOutput->Add(_h_eSFout_eSFin[j][i]);
	 fOutput->Add(_h_nPhe[j][i]);
      }
   }
   /* ****************************************************** */

   fOutput->Add(h_nParts);
   fOutput->Add(h_nPos);
   fOutput->Add(h_nNeg);
   fOutput->Add(h_nE);
   fOutput->Add(h_nP);
   fOutput->Add(h_nPip);
   fOutput->Add(h_nPim);
   fOutput->Add(h_nSurvive);
   fOutput->Add(h_nSurviveE_per);

}

Bool_t h10::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either h10::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.


   // fChain->GetEntry(entry);
   b_npart->GetEntry(entry);
   b_q_l->GetEntry(entry);
   b_stat->GetEntry(entry);
   b_p->GetEntry(entry);
   b_beta->GetEntry(entry);
   b_q->GetEntry(entry);
   b_id->GetEntry(entry);
   b_sc->GetEntry(entry);
   b_ec->GetEntry(entry);
   b_cc->GetEntry(entry);
   b_dc->GetEntry(entry);
   b_sc_sect->GetEntry(entry);
   b_ec_eo->GetEntry(entry);
   b_ec_ei->GetEntry(entry);
   b_nphe->GetEntry(entry);

   Int_t numPos = 0;
   Int_t numNeg = 0;
   Int_t numE = 0;
   Int_t numP = 0;
   Int_t numPip = 0;
   Int_t numPim = 0;
   Int_t charge, pid;

   for (int i = 0; i < npart; i++) {
      if (q_l == 0) return kFALSE;
      if (stat[i] < 0) continue;
      /* ******************************************************** */
      if (q[i] == -1 && sc[i] > 0 && ec[i] > 0 && cc[i] > 0 &&  dc[i] > 0) {
	 int scSect = sc_sect[sc[i]-1];
	 _h_eSF_p[0][scSect-1]->Fill(p[i], ec_eo[ec[i]-1]);
	 _h_eSFout_eSFin[0][scSect-1]->Fill(ec_ei[ec[i]-1], ec_eo[ec[i]-1]);
	 _h_nPhe[0][scSect-1]->Fill(nphe[cc[i]-1]);
	 if (id[i] == 11) {
	    _h_eSF_p[1][scSect-1]->Fill(p[i], ec_eo[ec[i]-1]);
	    _h_eSFout_eSFin[1][scSect-1]->Fill(ec_ei[ec[i]-1], ec_eo[ec[i]-1]);
	    _h_nPhe[1][scSect-1]->Fill(nphe[cc[i]-1]);
	 } else {
	    _h_eSF_p[3][scSect-1]->Fill(p[i], ec_eo[ec[i]-1]);
	    _h_eSFout_eSFin[3][scSect-1]->Fill(ec_ei[ec[i]-1], ec_eo[ec[i]-1]);
	    _h_nPhe[3][scSect-1]->Fill(nphe[cc[i]-1]);
	    if (id[i] == -211) {
	       _h_eSF_p[2][scSect-1]->Fill(p[i], ec_eo[ec[i]-1]);
	       _h_eSFout_eSFin[2][scSect-1]->Fill(ec_ei[ec[i]-1], ec_eo[ec[i]-1]);
	       _h_nPhe[2][scSect-1]->Fill(nphe[cc[i]-1]);
	    }
	 }
      }
      /* ******************************************************** */
      charge = q[i];
      pid = id[i];
      if (charge>0) numPos++;
      if (charge<0) {
	 numNeg++;
	 h_bvp[1]->Fill(p[i],beta[i]);
      }
      h_bvp[0]->Fill(p[i],beta[i]);
      int ispim = 0;
      switch(pid) {
      case 11:
	 numE++;
	 break;
      case 211:
	 numPip++;
	 break;
      case -211:
	 numPim++;
	 ispim = 1;
	 h_bvp[2]->Fill(p[i],beta[i]);
	 break;
      case 2212:
	 numP++;
	 break;
      }
      if (!ispim) h_bvp[3]->Fill(p[i],beta[i]);
   }

   h_nParts->Fill(npart);
   h_nPos->Fill(numPos);
   h_nNeg->Fill(numNeg);
   h_nE->Fill(numE);
   h_nP->Fill(numP);
   h_nPip->Fill(numPip);
   h_nPim->Fill(numPim);

   h_nSurvive->Fill(all,1);

   if ( numE>0 ) {
      h_nSurvive->Fill(n1e,1);
      nSurviveE++;
      nSurviveE_per++;
      if ( numP>0 ) {
	 h_nSurvive->Fill(n1e1p,1);
	 nSurviveEP++;
	 if ( numPip>0 ) {
	    h_nSurvive->Fill(n1e1p1pip,1);
	    nSurviveEPPip++;
	 }
	 if ( numPim>0 ) {
	    h_nSurvive->Fill(n1e1p1pim,1);
	    nSurviveEPPim++;
	 }
	 if ( numPip>0 && numPim>0 ) {
	    h_nSurvive->Fill(n1e1p1pip1pim,1);
	    nSurviveEPPipPim++;
	 }
      }
   }

   //ql_tmp = q_l;
   if (nfc == 0 || q_l != ql[nfc-1]) {
      ql[nfc] = q_l;
      if (nfc>0) {
	 //cout << "nfc = " << nfc << endl;
	 float qdiff = ql[nfc]-ql[nfc-1];
	 ql_sum += qdiff;
	 h_nSurviveE_per->Fill(nSurviveE_per/qdiff);
	 nSurviveE_per = 0;
      }
      nfc++;
   }
   if ( (numPos>=2 && numNeg>=1) || (numPos>=1 && numNeg>=2)) {
      if (numPos>=2 && numNeg>=1) h_nSurvive->Fill(n2pos1neg,1);
      if (numPos>=1 && numNeg>=2) h_nSurvive->Fill(n1pos2neg,1);
      if (numPos>=2 && numNeg>=2) h_nSurvive->Fill(n2pos2neg,1);
   }

   return kTRUE;
}

void h10::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.
   cout << "SlaveTerminate" << endl;
}

void h10::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
   cout << "Terminate" << endl;

   TFile *scratchFile = new TFile("scratch.root","recreate");
   scratchFile->cd();
   TIter next(fOutput);
   while (TObject *obj = next()) {
      if (obj->InheritsFrom("TH1")) {
	 //TH1 *tmpH = (TH1*)obj;
	 obj->Write(obj->GetName());
      }
   }

   //scratchFile->Close();
   const char *n_bvp[] = { "hbvp", "hbvp_neg", "hbvp_pim", "hbvp_pimNot" };
   TCanvas *c_bvp = new TCanvas("c_bvp","c_bvp");
   c_bvp->Divide(2,2);
   for (int i = 0; i < 4; i++) {
      c_bvp->cd(i+1);
      (fOutput->FindObject(n_bvp[i]))->Draw("colz");
   }

   TCanvas *c_eid[4];
   c_eid[0] = new TCanvas("c_eid1","c_eid");
   c_eid[1] = new TCanvas("c_eid2","c_eid only e^{-}");
   c_eid[2] = new TCanvas("c_eid3","c_eid only pi^{-}");
   c_eid[3] = new TCanvas("c_eid4","c_eid all except e^{-}");
   for (int i = 0; i < 4; i++) c_eid[i]->Divide(6,3);
   for(int i = 0; i < 6; i++){
      for(int j = 0; j < 4; j++) {
	 c_eid[j]->cd(i+1);
	 (fOutput->FindObject(_esf_p[j][i]))->Draw("colz");
	 gPad->SetLogz();
	 c_eid[j]->cd(6+i+1);
	 (fOutput->FindObject(_esf_inout[j][i]))->Draw("colz");
	 gPad->SetLogz();
	 c_eid[j]->cd(12+i+1);
	 (fOutput->FindObject(_nph[j][i]))->Draw();
      }
   }
}