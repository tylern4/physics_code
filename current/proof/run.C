
{  // set up a TChain
  TChain *ch = new TChain("h10", "My Chain");
  ch->Add("~/Data/e1d/new_cook/clas_*.root");

  plite = TProof::Open("");
  ch->SetProof();
  ch->Process("MySelector.cxx+");
  // gROOT->ProcessLine(".q;");
}
