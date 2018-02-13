
{  // set up a TChain
  TChain *ch = new TChain("h10", "My Chain");
  // r22853
  ch->Add("/Volumes/LaCiE/physics/e1d/v2/root/*.root");

  plite = TProof::Open("");
  ch->SetProof();
  ch->Process("MySelector.C+");
  gROOT->ProcessLine(".q;");
}
