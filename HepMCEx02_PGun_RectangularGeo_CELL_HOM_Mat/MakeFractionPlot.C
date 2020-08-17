

void MakeFractionPlot() {


  gStyle->SetOptStat(0);

  double EnergyVal[15] = {2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 30, 40, 50, 100};

  cout<<"Total = " << sizeof(EnergyVal)/sizeof(int) << endl;

  int NEnergy = (int)sizeof(EnergyVal)/sizeof(double);


  char name[500] ;

  double Efficiency[15], Resolution[15], InvSqrtE[15];

  for(int i_en = 0; i_en < 15; i_en++) {

  	sprintf(name, "PFlowNtupleFile_PiPlus%i.root", int(EnergyVal[i_en]) );
  	TFile *f = new TFile(name, "READ");

  	TTree *tr = (TTree*)f->Get("EventTree");
  	TH1D* ChEnergy = new TH1D("ChEnergy","ChEnergy",200,0,20000);
  	tr->Draw("Total_Ch_Energy>>ChEnergy", "", "HIST");

  	//TH1F *h = (TH1F*)f->Get("Total_Ch_Energy");
  	cout<<"Energy = " << EnergyVal[i_en] <<", Avg Mean = " << ChEnergy->GetMean()/(EnergyVal[i_en] * 100) << endl;

    Efficiency[i_en] = ChEnergy->GetMean()/(EnergyVal[i_en] * 100) ;
    Resolution[14 - i_en] = ChEnergy->GetRMS()/(EnergyVal[14 - i_en] * 100) ;
    InvSqrtE[14 - i_en] = 1/sqrt(EnergyVal[14 - i_en]);

  }

  TCanvas *c_resp = new TCanvas("HadronResponse", "HadronResponse", 600, 600);
  c_resp->cd();
  c_resp->SetLogx(1);
  c_resp->SetGridy(1);

  TGraph* gr = new TGraph(15, EnergyVal ,Efficiency);
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(1.4);
  gr->SetTitle("Response");
  gr->GetXaxis()->SetTitle("Energy [GeV]");
  gr->GetXaxis()->SetMoreLogLabels();
  gr->GetYaxis()->SetTitle("E_{reconstructed}/E_{available}");
  gr->GetYaxis()->SetRangeUser(0.5, 0.8);
  gr->Draw("AP");

  TCanvas *c_reso = new TCanvas("EnergyReso", "EnergyReso", 600, 600);
  c_reso->cd();
  //c_reso->SetLogx(1);
  c_reso->SetGridy(1);
 
  TGraph* gr_res = new TGraph(15, InvSqrtE ,Resolution);
  gr_res->SetMarkerStyle(20);
  gr_res->SetMarkerSize(1.4);
  gr_res->SetTitle("Response");
  gr_res->GetXaxis()->SetTitle("1/#sqrt{E}");
  gr_res->GetXaxis()->SetMoreLogLabels();
  gr_res->GetYaxis()->SetTitle("#frac{#sigma}{E}");
  //gr_res->GetYaxis()->SetRangeUser(0.5, 0.8);
  gr_res->Draw("AP");



 } 
