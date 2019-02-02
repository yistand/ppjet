void plotevent(int ievent = 350) {
	TFile *f = new TFile("NoTofMatch_FullJet_TransCharged_MatchTrig_ppJP_160811P12id_R06_HadrCorr_170127.root");
	TTree *ResultTree = (TTree*)f->Get("ResultTree");
	float jphi[100], jeta[100], sphi[100], seta[100], aphi[100], aeta[100], iphi[100], ieta[100];
	int j, s, a, i;
	float jp[100], sp[100], ap[100], ip[100];
	ResultTree->SetBranchAddress("TrkTranMinPt",ip);
	ResultTree->SetBranchAddress("TrkTranMaxPt",ap);
	ResultTree->SetBranchAddress("TrkSubAreaPt",sp);
	ResultTree->SetBranchAddress("TrkLeadAreaPt",jp);
	ResultTree->SetBranchAddress("LeadAreaNtrk",&j);
	ResultTree->SetBranchAddress("SubAreaNtrk",&s);
	ResultTree->SetBranchAddress("TranMaxNtrk",&a);
	ResultTree->SetBranchAddress("TranMinNtrk",&i);
	ResultTree->SetBranchAddress("TrkLeadAreaPhi",jphi);
	ResultTree->SetBranchAddress("TrkLeadAreaEta",jeta);
	ResultTree->SetBranchAddress("TrkSubAreaPhi",sphi);
	ResultTree->SetBranchAddress("TrkSubAreaEta",seta);
	ResultTree->SetBranchAddress("TrkTranMaxPhi",aphi);
	ResultTree->SetBranchAddress("TrkTranMaxEta",aeta);
	ResultTree->SetBranchAddress("TrkTranMinEta",ieta);
	ResultTree->SetBranchAddress("TrkTranMinPhi",iphi);

	ResultTree->GetEntry(ievent);

	gStyle->SetPalette(kDarkBodyRadiator);

	TH2D *h = new TH2D("h","",100,-3.14159,3.14159,20,-0.5,0.5);
	for(int k = 0; k<j; k++) {h->Fill(jphi[k],jeta[k],jp[k]);} for(int k = 0; k<s;k++) {h->Fill(sphi[k],seta[k],sp[k]);}  for(int k = 0; k<a; k++) {h->Fill(aphi[k],aeta[k],ap[k]);}  for(int k = 0; k<i; k++) {h->Fill(iphi[k],ieta[k],ip[k]);} 

	TCanvas *c1 = new TCanvas("c1","c1");
	h->Draw("lego2 psr");


	TCanvas *c2 = new TCanvas("c2","c2");
	TH1D *hpx = (TH1D*)h->ProjectionX("hpx");
	hpx->Draw("lego pol");







}
