#ifndef ROOT_DetaDphi
#define ROOT_DetaDphi


#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TString.h"
#include "TF1.h"

#include "Include/ClassTofMatchWeight.h"
#include "Include/ClassTPCWeight.h"

class DetaDphi{

	double MINPTCUT;	// CHECK!!!! need to cut consistence with production code in: 
				//UnderlyingAna::UnderlyingAna ( double R,
				//                double max_const_rap, //double PtConsLo, double PtConsHi,
				//                double min_const_pt,				<------------ This one
				//                double dPhiCut,
				//                TString name
				//                )
	double MINETACUT;	// |eta|<MINETACUT
	double EFFCUTOFF;	// if efficiency is too low, we will not use those particles to avoid too large correction

	double NeutralFracCut;

	int savefig;
	int saveroot;

	TFile *f;
	TTree *t;

	ClassTPCWeight *tpcweight;
	ClassTofMatchWeight *tofweight;

	TH2D *hetaphi;
	TH2D *hetaphi2;
	TH1D *hjet;
	TH1D *hjet2;

	TH1D *hjeteta;
	TH1D *hjet2eta;
	TH1D *hpeta;


	public:
	DetaDphi();
	~DetaDphi();
	
	float foldphi(float phi);
	float getweight(float eta, float pt, int charge, int MC);

	void MaxOrMin(float &max, float &min);
	void MaxOrMin(int &max, int &min);

	void SetSaveFig(int val);
	void SetSaveRoot(int val);

	void DoDetaDphi(TString dir, TString filetag, double jetptmin, double jetptmax, int ExclusiveEta, int DijetSameSide, double subjetptmin, double subjetptmax, int DoDijet);


};



#endif
