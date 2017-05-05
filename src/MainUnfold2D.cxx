//==============================================================================================================
//
//		2016.09.06	Li YI
//		use Unfold2D class to unfold undelrying vs leading jet pT
//
//
//==============================================================================================================

#if !(defined(__CINT__) || defined(__CLING__))
#include "Unfold2D.hh"

#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TFile.h>
#include <TString.h>
#include <TMath.h>
#include <TSystem.h>

#include <iostream>
#include <fstream>
#include <set>
#include <cmath>
#include <exception>
#include <cstdlib>      // std::rand, std::srand
#include <algorithm>    // std::random_shuffle

using namespace std;
#endif

void MainUnfold2D(TString ToUnfold="LeadAreaNtrk",const char* iter="5", const char* jetweight="1", const char* scaley = "0", const char *tpcsys="0", const char *tpcsyspm="0", const char *tpcsysabs="0", const char *Rc02Mc05="0") {
// tpcsys, tpcsyspm, tpcsysabs:
// if tpcsys==1, do TPC tracking efficiency uncertainty check: +/- (tpcsyspm) relative or absolute (tpcsysabs==0or1) 5%
// default tpcsys=0, don't apply TPC tracking sys. err. to unfolding training data
//
// Rc02Mc05:
// 	==1: Use RC pt>0.2, MC pt>0.5. Unfold to pt->0.5 using measured pt->0.2 data.
// 	==0: default, pt>0.2 for both RC and MC

	Unfold2D *uf2;


	//if((strcmp(iter,"4")!=0) || (strcmp(jetweight,"1")!=0) || (strcmp(scaley,"0")!=0) ) {
		const char *argv[16] = {"1","50","40","50","40","0","100","-0.5","39.5","0","1","1",iter,"1",jetweight, scaley};	// first opt: 1 for Bayes, 3 for Bin-by-Bin
		//  ntx, nty, nmx, nmy, xlo, xhi, ylo, yhi, overflow, verbose, doerror, regparm, WIDEBIN, flagjetweight, flagscaley
		uf2 = new Unfold2D("",16, argv);
	//}
	//else {
	//	uf2 = new Unfold2D("");
	//}
	
	uf2->Init();
	uf2->SetXYname("Leading jet p_{T}",ToUnfold);		// 2nd par: 	TranTotNtrk TranDiffNtrk TranTotPtSum TranPtAve NeutralFrac
									//		LeadAreaNtrk LeadPtAve LeadAreaPtSum
									//		SubAreaNtrk SubPtAve SubAreaPtSum
									//		JetNtrk, JetChargeNtrk
	uf2->SetTrigName("JP");						// MB JP0 JP1 JP2
	uf2->SetExcludeJPTrig(0);					// Embedding JP trigger exclude overlap or not. use 0 if unsure
	uf2->SetTranCharge("Charged");					// Charged Neutral for underlying particles
	//uf2->SetChangePrior(1);					// 1: for JP unfolding, use JP measured shape as prior for pythia MC level. 0: don't change pythia prior 
	uf2->SetTpcSys(atoi(tpcsys),atoi(tpcsyspm), atoi(tpcsysabs));			// TPC tracking 5% efficiency study, use absolute or relative +/- 5%

	uf2->SetRc02Mc05(atoi(Rc02Mc05));		// Use  measured pt>0.2 data to unfold back to pt>0.5

	uf2->PrintParms();
	//uf2->SetNoFake(true);
	//uf2->SetNoLoss(true);
/*
	uf2->TrainAndTest();
	uf2->Results();
	uf2->WriteTest();
*/



	//-- Use Train() TrainResults() WriteTrain() to read in MC data tree and fill histograms or use ReadResponseMatrix() if histogram is ready been fill by those functions into root file
	uf2->Train();
	uf2->TrainResults();
	uf2->WriteTrain();


/*
	// if already fill ResponseMatrix and root file available	
	uf2->ReadResponseMatrix();
*/
	


	//---- Use Fill4Unfold() to read in data tree and fill the histogram or use ReadMeasHist4Unfold() if the histogram is ready been fill by Fill4Unfold() in root file.
	// Recommend to use Fill4Unfold() instead of ReadMeasHist4Unfold(). THere is unsolved rebin to wide bin issue
	uf2->Fill4Unfold();
	//uf2->WriteHist4Unfold();		// if only want train histogram, don't do unfolding
	//uf2->ReadMeasHist4Unfold();
	uf2->Unfold();
	uf2->Results();
	uf2->WriteUnfoldResult();


	cout<<"Bye .. "<<endl;

}


#ifndef __CINT__
int main ( int argc, const char** argv) { 		// Main program when run stand-alone

	if( argc == 1) {
		MainUnfold2D(); 
	}
	else if( argc == 2) {
		vector<string> arguments(argv + 1, argv + argc);
		TString YName = arguments.at(0);
		MainUnfold2D(YName);
	}
	else if( argc == 3) {
		vector<string> arguments(argv + 1, argv + argc);
		TString YName = arguments.at(0);
		MainUnfold2D(YName, argv[2]);
	}
	else if( argc == 4) {
		vector<string> arguments(argv + 1, argv + argc);
		TString YName = arguments.at(0);
		MainUnfold2D(YName, argv[2], argv[3]);
	}
	else if( argc == 5) {
		vector<string> arguments(argv + 1, argv + argc);
		TString YName = arguments.at(0);
		MainUnfold2D(YName, argv[2], argv[3], argv[4]);
	}
	else if( argc == 6) {
		vector<string> arguments(argv + 1, argv + argc);
		TString YName = arguments.at(0);
		MainUnfold2D(YName, argv[2], argv[3], argv[4], argv[5]);
	}
	else if( argc == 7) {
		vector<string> arguments(argv + 1, argv + argc);
		TString YName = arguments.at(0);
		MainUnfold2D(YName, argv[2], argv[3], argv[4], argv[5], argv[6]);
	}
	else if( argc == 8) {
		vector<string> arguments(argv + 1, argv + argc);
		TString YName = arguments.at(0);
		MainUnfold2D(YName, argv[2], argv[3], argv[4], argv[5], argv[6], argv[7]);
	}
	else if( argc == 9) {
		vector<string> arguments(argv + 1, argv + argc);
		TString YName = arguments.at(0);
		MainUnfold2D(YName, argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8]);
	}

	return 0; 

}  
#endif

