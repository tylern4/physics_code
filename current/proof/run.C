
void run(std::string file = "~/Data/e1d/*.root") {  // set up a TChain
  TChain *ch = new TChain("h10", "My Chain");
  ch->Add(file.c_str());

  auto plite = TProof::Open("workers=16");
  ch->SetProof();
  ch->Process("MySelector.cxx");
  gROOT->ProcessLine(".q;");
}
