

void EnergyResponse() {

  gStyle->SetOptStat(0);

  TFile *f = TFile::Open("PFlowNtupleFile_QCD.root", "READ");
  
  TTree *tr = (TTree*)f->Get("EventTree");
  
  int NEntries = tr->GetEntries(); 

  cout << "Total Entries = " << NEntries << endl;


  TH1D* ChEnergy = new TH1D("ChEnergy","ChEnergy",100,0,2500);
  TH1D* NuEnergy = new TH1D("NuEnergy","NuEnergy",100,0,2500);


  tr->Draw("Total_Ch_Energy>>ChEnergy", "", "HIST");
  
  tr->Draw("Total_Nu_Energy>>NuEnergy", "", "HIST");

  TCanvas *canv = new TCanvas("EnergyDistribution", "EnergyDistribution", 600, 600);

  canv->cd();
  ChEnergy->GetXaxis()->SetTitle("Deposited Energy [GeV] X 100");
  ChEnergy->GetYaxis()->SetRangeUser(0, 1200);
  ChEnergy->SetLineWidth(3);
  ChEnergy->SetLineColor(2);
  ChEnergy->Draw("][");


  NuEnergy->SetLineWidth(3);
  NuEnergy->SetLineColor(4);
  NuEnergy->Draw("same ][");

 } // EnergyResponse()
