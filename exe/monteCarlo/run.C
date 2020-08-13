
void run(std::string file = "~/Data/sim/*.root") {  // set up a TChain
  TChain *ch = new TChain("h10", "My Chain");
  ch->Add(file.c_str());

  auto plite = TProof::Open("workers=4");
  ch->SetProof();
  ch->Process("sim.cxx+");
  // gROOT->ProcessLine(".q;");
}
