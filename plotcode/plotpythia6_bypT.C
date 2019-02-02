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

#include "CrossSectionPerpT.h"


void plotpythia6_bypT(bool flagMC05 = false) {
//
//	flagMC05:  whether we want to extract MC pt>0.5 from profile of MC pt>0.2 for TranPtAve
//

	const char *filenameformat = "pt%s_underMcVsEmbed_FullJetTransCharged_McPt02.root";

	//double ptbinning[13] = {0,2,3,4,5,7,9,11,15,20,25,35,45};
	const int nptbin = 13;		// need also to check outname..
	double ptbins[nptbin+1] = {0,1,3,5,7,9,11,15,20,25,35,45,55,100};		// 13+1 bins
	//const int nptbin = 100;
	//double ptbins[nptbin+1];
	//for(int i = 0; i<nptbin+1; i++) {
	//	ptbins[i] = i;
	//}


	//TH1D *hjpt = new TH1D("hjpt","leading jet pt",NUMBEROFPT,ptbinning);
	TH1D *hjpt = new TH1D("hjpt","leading jet pt",500,0,100);
	TProfile *pleadntrk = new TProfile("LeadNtrk","pleadntrk",nptbin,ptbins);
	TProfile *psubntrk = new TProfile("SubNtrk","psubntrk",nptbin,ptbins);
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
	pleadptsum->Sumw2();
	psubptsum->Sumw2();
	ptrantotptsum->Sumw2();
	pleadptave->Sumw2();
	psubptave->Sumw2();
	ptranptave->Sumw2();

	TH1D *htest = new TH1D("htest","htran_pt:8-9",20,0,20);

	int color[13] = {kOrange+10,kOrange+7,kOrange,kYellow+1,kSpring+9,kSpring-1,kGreen+3,kGreen-1,kTeal,kCyan+1,kAzure+3,kBlue-4};


	TList *listH = new TList();

	TList *listjpt = new TList();

	for(int i = 0; i<NUMBEROFPT; i++) {

		cout<<PTBINS[i]<<":"<<endl;

		TFile *f = new TFile(Form(filenameformat,PTBINS[i]));
		if(!f->IsOpen()) {cout<<"Cannot open "<<Form(filenameformat,PTBINS[i])<<endl; continue;}

		TTree *t = (TTree*)f->Get("ResultTree");
		float jpt;
		int leadntrk, subntrk, maxntrk, minntrk;
		float leadptsum, subptsum, maxptsum, minptsum;
		const int nmax=500;
		float leadpt[nmax],subpt[nmax],maxpt[nmax],minpt[nmax];
		double weight;

		t->SetBranchAddress("Mcj1pt",&jpt);
		//t->SetBranchAddress("eventweight",&weight);	// use the xsec from xsec file for each pT hat bin
        	t->SetBranchAddress("McLeadAreaNtrk",&leadntrk);
        	t->SetBranchAddress("McSubAreaNtrk",&subntrk);
        	t->SetBranchAddress("McTranMaxNtrk",&maxntrk);
        	t->SetBranchAddress("McTranMinNtrk",&minntrk);
        	t->SetBranchAddress("McLeadAreaPtSum",&leadptsum);
        	t->SetBranchAddress("McSubLeadAreaPtSum",&subptsum);
        	t->SetBranchAddress("McTranMaxPtSum",&maxptsum);
        	t->SetBranchAddress("McTranMinPtSum",&minptsum);
        	t->SetBranchAddress("McTrkLeadAreaPt",leadpt);
        	t->SetBranchAddress("McTrkSubAreaPt",subpt);
        	t->SetBranchAddress("McTrkTranMaxPt",maxpt);
        	t->SetBranchAddress("McTrkTranMinPt",minpt);

		TProfile *iptrantotntrk = new TProfile(Form("TranNtrk%s",PTBINS[i]),"iptrantotntrk",nptbin,ptbins);
		iptrantotntrk->Sumw2();

		TH1D *ihjpt = new TH1D(Form("ihjpt%s",PTBINS[i]),Form("leading jet pt %s",PTBINS[i]),500,0,100);


        	cout<<"Total # of Events: "<<t->GetEntries()<<endl;
        	for(int ievt = 0; ievt<t->GetEntries(); ievt++) {

			if(ievt%100000==0) cout<<ievt<<endl;		

        	        t->GetEntry(ievt);

			//weight = 1;	// TEST
			weight = XSEC[i]/NUMBEROFEVENT[i]; 		

			hjpt->Fill(jpt,weight);

			ihjpt->Fill(jpt,weight);

			if(!flagMC05) {
				pleadntrk->Fill(jpt,leadntrk,weight);		
				psubntrk->Fill(jpt,subntrk,weight);		
				ptrantotntrk->Fill(jpt,maxntrk+minntrk,weight);		

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

				iptrantotntrk->Fill(jpt,maxntrk05+minntrk05,weight);		

				pleadptsum->Fill(jpt,leadptsum05,weight);		
				psubptsum->Fill(jpt,subptsum05,weight);		
				ptrantotptsum->Fill(jpt,maxptsum05+minptsum05,weight);		

				if(jpt<9 && jpt>8) htest->Fill(maxntrk05+minntrk05,weight);
			
			}
		}
		cout<<PTBINS[i]<<": "<<iptrantotntrk->GetBinContent(5)<<" "<<iptrantotntrk->GetBinEntries(5)<<" "<<iptrantotntrk->GetBinEffectiveEntries(5)<<endl;
		listH->Add(iptrantotntrk);
		listjpt->Add(ihjpt);
	}

	//TString outname = "Profile"+filename;
	//TString outname = "Profile_finebin_pythia6.root";
	TString outname = "Profile_12JetBinv2_pythia6.root";
	if(flagMC05) outname.ReplaceAll(".root","_PT05.root");
	//outname.ReplaceAll(".root","_noxsecw.root");	// TEST
	TFile *fout = new TFile(outname,"RECREATE");
	hjpt->Write();
	pleadntrk->Write();
	psubntrk->Write();
	ptrantotntrk->Write();
	pleadptsum->Write();
	psubptsum->Write();
	ptrantotptsum->Write();
	pleadptave->Write();
	psubptave->Write();
	ptranptave->Write();
	htest->Write();
	listH->Write();
	listjpt->Write();

	fout->Close();

/*
	TCanvas *c = new TCanvas();
	gStyle->SetOptStat(0);
	c->SetLogy();

	hjpt->GetXaxis()->SetTitle("Leading jet p_{T}");
	hjpt->GetYaxis()->SetTitle("#sigma");
	hjpt->Draw();

	TIter next(listjpt);
	TObject *obj;
	int icount = 0;
	while((obj = next.Next())) {
		TH1D *htmp = dynamic_cast<TH1D*>(obj);
		//cout<<htmp<<endl;
		if(!htmp) continue;
		htmp->SetLineColor(color[icount]);
		htmp->SetMarkerColor(color[icount]);
		icount++;
		htmp->Draw("HISTsame");
	}
*/

}
