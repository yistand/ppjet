//===============================================================================================================
//
//	2016.03.22	Li Yi
//	not sure why ppjet code I run today (I don't recall any major change will cause it). 
//
//
//===============================================================================================================

#include <iostream>
#include <vector>
#include <algorithm>

#include "TFile.h"
#include "TH2D.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TProfile.h"
#include "TMath.h"

#include "ResultTree.C"
#include "/lustre/home/client/hep/caines/ly247/ppjet/plotcode/Jet.C"

using namespace std;



vector<int> GetEventid(int runid, const char *inputvector="FullJet_TransCharged_MatchTrig_ppJP2_151030P12id_R06_HadrCorr_160322.root") {

        TChain *chain = new TChain("ResultTree");
	chain->Add(inputvector);

	ResultTree *t = new ResultTree(chain);

	vector<int> eventid;
	
	for(int ievt = 0 ; ievt<chain->GetEntries() ; ievt++) {
		t->GetEntry(ievt);
		if(t->runid==runid) {
			eventid.push_back(t->eventid);
		}
	}

	sort(eventid.begin(), eventid.end());
	cout<<"event size = "<<eventid.size();
	unique(eventid.begin(),eventid.end());
	cout<<" after deleting duplicate elements -> "<<eventid.size()<<endl;

	return eventid;

}


vector<int> compare2(int runid=13071005, const char *input1="FullJet_TransCharged_MatchTrig_ppJP2_151030P12id_R06_HadrCorr_160314.root", const char *input2="FullJet_TransCharged_MatchTrig_ppJP2_151030P12id_R06_HadrCorr_160322.root")	{
	
	vector<int> v314 = GetEventid(runid,input1);
	vector<int> v322 = GetEventid(runid,input2);

	for(unsigned int i = 0; i<v322.size(); i++) {
		for(unsigned int j= 0; j<v314.size(); j++) {
			if(v314.at(j)==v322.at(i)) v314.erase(v314.begin()+j);
		}
	}

	cout<<"delete v322 from v314: "<<v314.size()<<endl;

	vector<int>::iterator it;
	for(it=v314.begin(); it!=v314.end(); ++it) {
		cout<<" "<<*it;
	}
	cout<<endl;

	return v314;

}


void plotVpdVz(int runid=13071005, const char *jetinput = "/home/hep/caines/ly247/Scratch/pp12JP2Pico_151030/*2.root") { 

	TH1D *h = new TH1D("h","missed vpdVz",2000,-1000,1000);
	h->SetMarkerStyle(8);
	TH1D *h2 = new TH1D("h2","passed vpdVz",2000,-1000,1000);
	h2->SetMarkerStyle(8);

	TChain *chain = new TChain("JetTree");
	chain->Add(jetinput);

	Jet *t = new Jet(chain);

	vector<int> missing = compare2(runid);
	vector<int> passed = GetEventid(runid);
	vector<int>::iterator it;
	for(int ievt = 0 ; ievt<chain->GetEntries() ; ievt++) {
		t->GetEntry(ievt);
		if(t->fEventHeader_fRunId==runid) {
			it = find(missing.begin(),missing.end(),t->fEventHeader_fEventId);
			if(it!=missing.end())	{
				cout<<(t->fEventHeader_fRunId)<<" "<<(t->fEventHeader_fEventId)<<" "<<(t->fEventHeader_fvpdVz)<<" diff: "<<(t->fEventHeader_fvpdVz)-(t->fEventHeader_fPVz)<<endl;
				h->Fill(t->fEventHeader_fvpdVz);
			}
			else {
				it = find(passed.begin(),passed.end(),t->fEventHeader_fEventId);
				if(it!=passed.end()) {
					h2->Fill(t->fEventHeader_fvpdVz);
				}
			}
		}
	}

	TCanvas *c1 = new TCanvas();
	h->Draw();

	TCanvas *c2 = new TCanvas();
	h2->Draw();



	
}

