#define MIN(a, b) ((a < b) ? a : b)
#define MAX(a, b) ((a > b) ? a : b)

double poissonf(double *x, double *par) { return par[0] * TMath::Poisson(x[0], par[1]); }

void run() {
  const char *fin = "5ms_1000s.root";

  TChain *chain = new TChain("Sr90");
  chain->Add(fin);
  int _x1;
  chain->SetBranchAddress("x1", &_x1);

  int max = 0;
  int min = 100000;
  int num_of_events = (int)chain->GetEntries();
  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain->GetEntry(current_event);
    max = MAX(max, _x1);
    min = MIN(min, _x1);
  }

  int bins = num_of_events / 1000;
  max = max + max * 0.1;
  min = min - min * 0.1;
  std::cout << "Making histogram with " << bins;
  std::cout << " bins and max/min values of (";
  std::cout << max << "," << min << ")" << std::endl;
  TH1I *hist = new TH1I("x1", "x1", bins, min, max);

  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain->GetEntry(current_event);
    hist->Fill(_x1);
  }

  TF1 *pois = new TF1("pois", poissonf, 0, 10, 2);

  pois->SetParName(0, "Const");
  pois->SetParName(1, "#mu");
  pois->SetParameter(0, 1);
  pois->SetParameter(1, 1);

  hist->Fit("pois");
  chain->Reset();
  hist->Draw();
  gStyle->SetOptStat(111);
  gStyle->SetOptFit(111);
}
