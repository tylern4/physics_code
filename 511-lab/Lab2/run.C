#define MIN(a, b) ((a < b) ? a : b)
#define MAX(a, b) ((a > b) ? a : b)

double poissonf(double *x, double *par) { return par[0] * TMath::Poisson(x[0], par[1]); }
double gausf(double *x, double *par) { return par[0] * TMath::Gaus(x[0], par[1], par[2], par[3]); }
double line(double *x, double *par) { return par[0] + x[0] * par[1]; }

const char *fin[8] = {"0.5ms_1000s", "0.5ms_100s", "0.5ms_10s", "50ms_10000s",
                      "5ms_1000s",   "5ms_100s",   "5ms_10s",   "500ms_50000s.n0"};

TCanvas c1("c1", "Plots", 1280, 720);
TCanvas c2("c2", "Graph", 640, 360);

const int bins[8] = {8, 8, 8, 65, 17, 16, 13, 125};
const int min_val[8] = {0, 0, 0, 25, 0, 0, 0, 425};
const int max_val[8] = {7, 5, 5, 90, 17, 16, 13, 675};
const float dt[8] = {0.5, 0.5, 0.5, 50, 5, 5, 5, 500};
const float T[8] = {1000, 100, 10, 10000, 1000, 100, 10, 50000};

float process(int file) {
  // ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  TChain *chain = new TChain("Sr90");
  char file_name[50];
  sprintf(file_name, "%s.root", fin[file]);
  chain->Add(file_name);

  int _x1;
  chain->SetBranchAddress("x1", &_x1);
  int num_of_events = (int)chain->GetEntries();
  TH1I *hist = new TH1I(fin[file], fin[file], bins[file], min_val[file], max_val[file]);
  for (int current_event = 0; current_event < num_of_events; current_event++) {
    chain->GetEntry(current_event);
    hist->Fill(_x1);
  }

  chain->Reset();
  c1.cd((int)file + 1);
  hist->SetXTitle("x");
  hist->Draw();

  TF1 *pois = new TF1("pois", poissonf, min_val[file], max_val[file], 2);
  pois->SetParName(0, "Const");
  pois->SetParName(1, "#mu");
  pois->SetParameter(0, 1);
  pois->SetParameter(1, hist->GetMean());
  hist->Fit("pois", "QM+", "", min_val[file], max_val[file]);
  return pois->GetParameter("#mu");
}

void graph(float mu[8]) {
  TGraph *gr = new TGraph(8, mu, dt);
  TF1 *lines = new TF1("lines", line, 0, 1000, 2);
  lines->SetParName(0, "Intercept");
  lines->SetParName(1, "Slope");
  c2.cd();
  gr->Fit("lines", "QM+", "", 0, 1000);
  gr->Draw("A*");
}

void run() {
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(1111);
  std::vector<float> r_vec;
  float mu[8];
  c1.Divide(4, 2);
  for (int i = 0; i < 8; i++) {
    mu[i] = process(i);
  }
  std::cout << "mu,dt,T,r" << std::endl;
  for (int i = 0; i < 8; i++) std::cout << mu[i] << "," << dt[i] << "," << T[i] << "," << mu[i] / dt[i] << std::endl;
  gStyle->SetStatX(0.5);
  gStyle->SetStatY(0.8);
  graph(mu);
}
