//=====================================================================================================
//
//	2017.11.15	Li YI
//	get pythia8 leading jet pt by scale each pt bin with cross section/nEvents
//
//=====================================================================================================


void PlotPythia8LeadJetpT() {

	const char *tag = "171114";

	const char *tag2 = ""; //"_MPIoff";

	TString inputdir = "./"; //"/home/fas/caines/ly247/Scratch/pythiadata/";
	const char* inputfiletag = "FullJet_TransCharged_pythia8215_pp200hard%s_pT%d-%d_PionDecayOff_seed%s_170422.root";

	const char* xsecfile = "Xsec4pythia8215_pp200hard%s_PionDecayOff_seed%s.csv";

	int color[13] = {kRed,kOrange+10,kOrange+7,kOrange,kYellow+1,kSpring+9,kSpring-1,kGreen+3,kGreen-1,kTeal,kCyan+1,kAzure+3,kBlue-4};

	ifstream fin;
	fin.open(inputdir+Form(xsecfile,tag2,tag));
	if(!fin.is_open()) {cout<<"Cannot open "<<Form(xsecfile,tag2,tag)<<endl; return;}

	string line;
	float pt1 = 0, pt2 = 0;
	float xsec = 0, evts = 0;

	TH1D *hsum = new TH1D("hsum","Leading jet p_{T}",500,0,100);
	TList *listH = new TList();

	while(getline(fin, line)) {		// loop each pT bin
		//cout<<line<<endl;
		if ( line.size()==0 ) continue; // skip empty lines
		if ( line[0] == '#' ) continue; // skip comments

		std::istringstream ss( line );

		char ctmp;

		ss >> pt1 >> ctmp >> pt2 >> ctmp >> xsec >> ctmp >> evts;
		cout<<pt1<<" "<<pt2<<" "<<xsec<<" "<<evts<<endl;

		TFile *f1 = new TFile(inputdir+Form(inputfiletag,tag2,int(pt1),int(pt2),tag));
		if(!f1->IsOpen()) {cout<<"Cannot open "<<Form(inputfiletag,tag2,int(pt1),int(pt2),tag)<<endl; return;}

		float jpt;
		TTree *t = (TTree*)f1->Get("ResultTree");
		if(!t) {cout<<"Cannot find TTree t in "<<Form(inputfiletag,int(pt1),int(pt2),tag)<<endl; return;}
		t->SetBranchAddress("j1pt",&jpt);

		TH1D *hpti = new TH1D(Form("hpt%d_%d",int(pt1),int(pt2)),Form("Leading jet p_{T} for p_{T}^{Hat}: %d-%d",int(pt1),int(pt2)),500,0,100);

		for(int ievt = 0; ievt<t->GetEntries(); ievt++) {

			t->GetEntry(ievt);

			hpti->Fill(jpt);
		}

		hpti->Scale(xsec/evts);
		hsum->Add(hpti);

		listH->Add(hpti);

	}

	TCanvas *c = new TCanvas();
	gStyle->SetOptStat(0);
	c->SetLogy();

	hsum->GetXaxis()->SetTitle("Leading jet p_{T}");
	hsum->GetYaxis()->SetTitle("#sigma");
	hsum->Draw();

	hsum->SetTitle(tag2);

	listH->ls();
	TIter next(listH);
	TObject *obj;
	int icount = 0;
	while((obj = next.Next())) {
		TH1D *htmp = dynamic_cast<TH1D*>(obj);
		htmp->SetLineColor(color[icount]);
		htmp->SetMarkerColor(color[icount]);
		icount++;
		htmp->Draw("HISTsame");
	}

}
