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

void MainUnfold2D(TString ToUnfold="JetChargeNtrk",const char* iter="4", const char* jetweight="1", const char* scaley = "0") {

	Unfold2D *uf2;


	if((strcmp(iter,"4")!=0) || (strcmp(jetweight,"1")!=0) || (strcmp(scaley,"0")!=0) ) {
		const char *argv[16] = {"1","50","40","50","40","0","100","-0.5","39.5","0","1","1",iter,"1",jetweight, scaley};
		uf2 = new Unfold2D("",16, argv);
	}
	else {
		uf2 = new Unfold2D("");
	}
	
	uf2->Init();
	uf2->SetXYname("Leading jet p_{T}",ToUnfold);		// 2nd par: 	TranTotNtrk TranTotPtSum TranPtAve NeutralFrac
									//		LeadAreaNtrk
									//		SubAreaNtrk
									//		JetNtrk, JetChargeNtrk
	uf2->SetTrigName("JP");						// MB JP0 JP1 JP2
	uf2->SetExcludeJPTrig(0);					// Embedding JP trigger exclude overlap or not
	uf2->SetTranCharge("Charged");					// Charged Neutral

	uf2->PrintParms();
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

	return 0; 

}  
#endif

