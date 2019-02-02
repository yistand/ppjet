//===================================================================================================================
//
//	2017.05.10	Li YI
//	Read ResultTree from pythia 8 file FullJet_TransCharged_pythia8215_pp200hard_PionDecayOff_seed134123_170422.root
//	(decay off for pid 111,211,221,321,310,130,3122,3212,3112,3222,3322,3334)
//	Output: 
//	Lead, Sub, Tran,   PtAve, Ntrk, PtSum, ProfileX() vs lead jet pt
//
//===================================================================================================================

#include "TH1D.h"
#include "TProfile.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"




//void plotpythia8_bypT(const char *tag = "", const char *var = "maxpt", bool flagMC05 = false) {	
void plotpythia8_bypT(const char *tag = "", const char *var = "j1pt", bool flagMC05 = false) {	
//
//	flagMC05:  whether we want to extract MC pt>0.5 from profile of MC pt>0.2 for TranPtAve
//
//	var: maxpt for MaxTrack, j1pt for Leading jet
//
//	tag: "" or "_MPIoff"

	const char *filenameformat = "%sFullJet_TransCharged_pythia8215_pp200hard%s_pT%d-%d_PionDecayOff_seed171114_170422.root";
	//const char *filenameformat = "%sChargeJet_TransCharged_pythia8215_pp200hard%s_pT%d-%d_PionDecayOff_seed171114_170422.root";
	//const char *filenameformat = "%sCMSsetting_FullJet_TransCharged_pythia8215_pp200hard%s_pT%d-%d_PionDecayOff_seed171114_170422.root";

	const char* xsecfile = "Xsec4pythia8215_pp200hard%s_PionDecayOff_seed171114.csv";


	ifstream xfin;
	xfin.open(Form(xsecfile,tag));
	if(!xfin.is_open()) {cout<<"Cannot open "<<Form(xsecfile,tag)<<endl; return;}

	string line;
	float pt1 = 0, pt2 = 0;
	float xsec = 0, evts = 0;

	//const int nptbin = 20;		// need also to check outname..
	//double ptbins[nptbin+1] = {0,2,3,4,5,6,7,8,9,11,12,15,17,20,22,25,30,35,45,65,100};	// 20+1 bins
	const int nptbin = 13;		// need also to check outname..
	double ptbins[nptbin+1] = {0,1,3,5,7,9,11,15,20,25,35,45,55,100};		// 13+1 bins
	//const int nptbin = 100;
	//double ptbins[nptbin+1];
	//for(int i = 0; i<nptbin+1; i++) {
	//	ptbins[i] = i;
	//}

	TH1D *hjpt = new TH1D("hjpt","leading jet pt",500,0,100);
	TProfile *pleadntrk = new TProfile("LeadNtrk","pleadntrk",nptbin,ptbins);
	TProfile *psubntrk = new TProfile("SubNtrk","psubntrk",nptbin,ptbins);
	TProfile *ptranmaxntrk = new TProfile("TranMaxNtrk","ptranmaxntrk",nptbin,ptbins);
	TProfile *ptranminntrk = new TProfile("TranMinNtrk","ptranminntrk",nptbin,ptbins);
	TProfile *ptrantotntrk = new TProfile("TranNtrk","ptrantotntrk",nptbin,ptbins);
	TProfile *pleadptsum = new TProfile("LeadPtSum","pleadptsum",nptbin,ptbins);
	TProfile *psubptsum = new TProfile("SubPtSum","psubptsum",nptbin,ptbins);
	TProfile *ptrantotptsum = new TProfile("TranPtSum","ptrantotptsum",nptbin,ptbins);
	TProfile *pleadptave = new TProfile("LeadPtAve","pleadptave",nptbin,ptbins);
	TProfile *psubptave = new TProfile("SubPtAve","psubptave",nptbin,ptbins);
	TProfile *ptranptave = new TProfile("TranPtAve","ptranptave",nptbin,ptbins);

	pleadntrk->Sumw2();
	psubntrk->Sumw2();
	ptrantotntrk->Sumw2();
	ptranmaxntrk->Sumw2();
	ptranminntrk->Sumw2();
	pleadptsum->Sumw2();
	psubptsum->Sumw2();
	ptrantotptsum->Sumw2();
	pleadptave->Sumw2();
	psubptave->Sumw2();
	ptranptave->Sumw2();

	TH1D *htest = new TH1D("htest","htran_pt:8-9",20,0,20);

	TList *listH = new TList();

	TString pretag = "";
	if(!strcmp(var,"maxpt")) pretag = "MaxTrack_";
	while(getline(xfin, line)) {		// loop each pT bin
		//cout<<line<<endl;
		if ( line.size()==0 ) continue; // skip empty lines
		if ( line[0] == '#' ) continue; // skip comments

		std::istringstream ss( line );

		char ctmp;

		ss >> pt1 >> ctmp >> pt2 >> ctmp >> xsec >> ctmp >> evts;
		cout<<pt1<<" "<<pt2<<" "<<xsec<<" "<<evts<<endl;

		TFile *f = new TFile(Form(filenameformat, pretag.Data(), tag, int(pt1),int(pt2)));
		if(!f->IsOpen()) {cout<<"Cannot open "<<Form(filenameformat, pretag.Data(), tag, int(pt1),int(pt2))<<endl; continue;}
		else {cout<<"Read in "<<f->GetName()<<endl;}

		TTree *t = (TTree*)f->Get("ResultTree");
		float jpt;
		int leadntrk, subntrk, maxntrk, minntrk;
		float leadptsum, subptsum, maxptsum, minptsum;
		const int nmax=500;
		float leadpt[nmax],subpt[nmax],maxpt[nmax],minpt[nmax];
		double weight;

		t->SetBranchAddress(var,&jpt);		// "maxpt" for track, or "j1pt" for leading jet
		//t->SetBranchAddress("j1pt",&jpt);
		//t->SetBranchAddress("eventweight",&weight);	// use the xsec from xsec file for each pT hat bin
        	t->SetBranchAddress("LeadAreaNtrk",&leadntrk);
        	t->SetBranchAddress("SubAreaNtrk",&subntrk);
        	t->SetBranchAddress("TranMaxNtrk",&maxntrk);
        	t->SetBranchAddress("TranMinNtrk",&minntrk);
        	t->SetBranchAddress("LeadAreaPtSum",&leadptsum);
        	t->SetBranchAddress("SubLeadAreaPtSum",&subptsum);
        	t->SetBranchAddress("TranMaxPtSum",&maxptsum);
        	t->SetBranchAddress("TranMinPtSum",&minptsum);
        	t->SetBranchAddress("TrkLeadAreaPt",leadpt);
        	t->SetBranchAddress("TrkSubAreaPt",subpt);
        	t->SetBranchAddress("TrkTranMaxPt",maxpt);
        	t->SetBranchAddress("TrkTranMinPt",minpt);

		TProfile *iptrantotntrk = new TProfile(Form("TranNtrk%d_%d",int(pt1),int(pt2)),"iptrantotntrk",nptbin,ptbins);

		iptrantotntrk->Sumw2();

        	cout<<"Total # of Events: "<<t->GetEntries()<<endl;
        	for(int ievt = 0; ievt<t->GetEntries(); ievt++) {

			if(ievt%100000==0) cout<<ievt<<endl;		

        	        t->GetEntry(ievt);

			weight = xsec/evts; 			
			//weight = 1; 			// TEST

			hjpt->Fill(jpt,weight);

			if(!flagMC05) {
				pleadntrk->Fill(jpt,leadntrk,weight);		
				psubntrk->Fill(jpt,subntrk,weight);		
				ptrantotntrk->Fill(jpt,maxntrk+minntrk,weight);		
				ptranmaxntrk->Fill(jpt,maxntrk,weight);		
				ptranminntrk->Fill(jpt,minntrk,weight);		

				iptrantotntrk->Fill(jpt,maxntrk+minntrk,weight);		

				pleadptsum->Fill(jpt,leadptsum,weight);		
				psubptsum->Fill(jpt,subptsum,weight);		
				ptrantotptsum->Fill(jpt,maxptsum+minptsum,weight);		
			
				for(int i = 0; i<leadntrk; i++) {
					pleadptave->Fill(jpt,leadpt[i],weight);
				}
				for(int i = 0; i<subntrk; i++) {
					psubptave->Fill(jpt,subpt[i],weight);
				}
				for(int i = 0; i<maxntrk; i++) {
					ptranptave->Fill(jpt,maxpt[i],weight);
				}
				for(int i = 0; i<minntrk; i++) {
					ptranptave->Fill(jpt,minpt[i],weight);
				}

				if(jpt<9 && jpt>8) htest->Fill(maxntrk+minntrk,weight);
			}
			else {
				int leadntrk05 = 0, subntrk05 = 0, maxntrk05 = 0, minntrk05 = 0;
				float leadptsum05 = 0, subptsum05 = 0, maxptsum05 = 0, minptsum05 = 0;
				for(int i = 0; i<leadntrk; i++) {
					if(leadpt[i]>0.5) {
						pleadptave->Fill(jpt,leadpt[i],weight);
						leadntrk05++;
						leadptsum05+=leadpt[i];
					}
				}
				for(int i = 0; i<subntrk; i++) {
					if(subpt[i]>0.5) {
						psubptave->Fill(jpt,subpt[i],weight);
						subntrk05++;
						subptsum05+=subpt[i];
					}
				}
				for(int i = 0; i<maxntrk; i++) {
					if(maxpt[i]>0.5) {
						ptranptave->Fill(jpt,maxpt[i],weight);
						maxntrk05++;
						maxptsum05+=maxpt[i];
					}
		
				}
				for(int i = 0; i<minntrk; i++) {
					if(minpt[i]>0.5) {
						ptranptave->Fill(jpt,minpt[i],weight);
						minntrk05++;
						minptsum05+=minpt[i];
					}
				}
				pleadntrk->Fill(jpt,leadntrk05,weight);		
				psubntrk->Fill(jpt,subntrk05,weight);		
				ptrantotntrk->Fill(jpt,maxntrk05+minntrk05,weight);		
				ptranmaxntrk->Fill(jpt,maxntrk05,weight);		
				ptranminntrk->Fill(jpt,minntrk05,weight);		

				iptrantotntrk->Fill(jpt,maxntrk05+minntrk05,weight);		

				pleadptsum->Fill(jpt,leadptsum05,weight);		
				psubptsum->Fill(jpt,subptsum05,weight);		
				ptrantotptsum->Fill(jpt,maxptsum05+minptsum05,weight);		

				if(jpt<9 && jpt>8) htest->Fill(maxntrk05+minntrk05,weight);
			
			}
		}
		cout<<pt1<<"-"<<pt2<<": "<<iptrantotntrk->GetBinContent(5)<<" "<<iptrantotntrk->GetBinEntries(5)<<" "<<iptrantotntrk->GetBinEffectiveEntries(5)<<endl;

		listH->Add(iptrantotntrk);
	}

	//TString outname = pretag+Form("Profile_finebin_pythia8%s_pp200hard_PionDecayOff_seed171114_170422.root",tag);;
	TString outname = pretag+Form("Profile_12JetBinv2_pythia8%s.root",tag);
	if(flagMC05) outname.ReplaceAll(".root","_PT05.root");
	if(strstr(filenameformat,"ChargeJet")!=NULL) outname.ReplaceAll("Profile","ProfileVsChargeJet");
	//outname.ReplaceAll(".root","_noxsecw.root");	// TEST
	TFile *fout = new TFile(outname,"RECREATE");
	hjpt->Write();
	pleadntrk->Write();
	psubntrk->Write();
	ptrantotntrk->Write();
	ptranmaxntrk->Write();
	ptranminntrk->Write();
	pleadptsum->Write();
	psubptsum->Write();
	ptrantotptsum->Write();
	pleadptave->Write();
	psubptave->Write();
	ptranptave->Write();
	htest->Write();
	listH->Write();

	fout->Close();

}
