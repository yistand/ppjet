#ifndef ROOT_RCVZWEIGHT
#define ROOT_RCVZWEIGHT

#include "TH1D.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"

class RcVzWeight{
	
	TF1 *fpol;	// fit function

	TH1D *hw;	// ratio histogram

	double VzRange;	// -30 -> 30 cm
	
	TH1D *hrc;
	TH1D *hdata;

	TFile *frc;
	TTree *treerc;

	int runidrc;
	float jptrc;
	double vzrc;

	double xsecw;

	TFile *fdata;
	TTree *treedata;

	int runiddata;
	float jptdata;
	double vzdata;


	public:
	RcVzWeight();
	~RcVzWeight();

	double GetVzWeight(double vz);
	float GetVzWeight(float vz);
	bool ReadAndFill(TString filetagrc="JP", TString filenamedata="/home/fas/caines/ly247/Scratch/pp200Y12_jetunderlying/NoTofMatch_FullJet_TransCharged_MatchTrig_ppJP_160811P12id_R06_HadrCorr_170418.root");
	bool WriteVzFile(TString outfilename="VzWeight.root");
	bool ReadVzfile(TString infilename="VzWeight.root");
	bool SetVzRange(double range);
	bool InitHist(TString filetagrc);
	

};

#endif



