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

void MainUnfold2D() {

	Unfold2D *uf2 = new Unfold2D("trig");
	
	uf2->Init();
	uf2->SetXYname("Leading jet p_{T}","TranTotNtrk");
	uf2->SetTrigName("JP1");

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


	//uf2->ReadResponseMatrix();
	
	//---- Use Fill4Unfold() to read in data tree and fill the histogram or use ReadMeasHist4Unfold() if the histogram is ready been fill by Fill4Unfold() in root file.
	//uf2->Fill4Unfold();
	uf2->ReadMeasHist4Unfold();
	uf2->Unfold();
	uf2->Results();
	uf2->WriteUnfoldResult();

	cout<<"Bye .. "<<endl;

}


#ifndef __CINT__
int main () { MainUnfold2D(); return 0; }  // Main program when run stand-alone
#endif

