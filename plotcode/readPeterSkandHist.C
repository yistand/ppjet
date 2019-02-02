void readPeterSkandHist() {

	TGraphErrors *gr = new TGraphErrors();

	ifstream fin;
	fin.open("TranNtrkVsJptPeterSkandpp7TeVMB.dat");
	float pt1, pt2;
	float y, ey1, ey2;
	string line;
	int ict;

	while(getline(fin,line)) {
		istringstream ss(line);
		ss >> ict >> pt1 >> pt2 >> y >> ey1 >> ey2;
		ict=ict-1;
		gr->SetPoint(ict, pt1, y);
		gr->SetPointError(ict, 0, ey1);
		cout<<pt1<<" "<<y<<" "<<ey1<<endl;
	}

	TFile *f = new TFile("R0.5_ChargeJet_TransCharged_pythia8215_pp7TeVMB_seed105892_171120.root");
	TProfile *ptran = (TProfile*)f->Get("TranNtrkvsLeadJetPt");
	float DeDpNorma = 1./(4.*TMath::Pi()/3.);	// TranAve: eta 4; phi pi/3
	ptran->Scale(DeDpNorma);

	TFile *f_star = new TFile("FullJet_TransCharged_pythia8215_pp7TeVMB_seed105892_170422_STARsetting.root");
	TProfile *ptran_star = (TProfile*)f_star->Get("TranNtrkvsLeadJetPt");
	float DeDpNorma_star = 1./(2.*TMath::Pi()/3.);	// TranAve: eta 2; phi pi/3
	ptran_star->Scale(DeDpNorma_star);
	ptran_star->Scale(0.5);

	ptran->GetXaxis()->SetTitle("Leading Charged #color[3]{(Full)} Jet p_{T}");
	ptran->GetYaxis()->SetTitle("Transverse N/d#etad#phi");
	ptran->SetTitle("Transverse Charged Density vs Leading Charged #color[3]{(Full)} Jet p_{T}");
	ptran->GetXaxis()->SetRangeUser(0,45);
	ptran->GetYaxis()->SetRangeUser(0,2);

	ptran->SetMarkerStyle(25);
	ptran->SetMarkerColor(1);
	ptran->SetLineColor(1);

	gr->SetMarkerStyle(4);
	gr->SetMarkerColor(2);
	gr->SetLineColor(2);

	ptran_star->SetMarkerStyle(34);
	ptran_star->SetMarkerSize(1.5);
	ptran_star->SetMarkerColor(kGreen);
	ptran_star->SetLineColor(kGreen);

	ptran->Draw();
	ptran_star->Draw("epsame");
	gr->Draw("epsame");

	gStyle->SetOptStat(0);

	TLegend *leg = new TLegend(0.15,0.7,0.6,0.85);
	leg->AddEntry(gr,"Peter's","pl");
	leg->AddEntry(ptran,"Li's","pl");
	leg->AddEntry(ptran_star,"Li's STAR setting (x0.5)","pl");
	leg->Draw("same");

	TLatex *lat = new TLatex(4,1.2,"p+p@7 TeV");
	lat->SetTextSize(0.08);
	lat->Draw("same");
}
