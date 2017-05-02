#ifndef ROOT_REWEIGHTBYY
#define ROOT_REWEIGHTBYY

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TProfile.h"

class ScaleByY{

	TString YVariableName;
	
	TH1F *h1w;

	TProfile *pfjp;
	TProfile *pfmb;


	TH1D *htmprx;

	TFile *fjp;
	TTree *treejp;

	float Yjp;
	float Ymb;

	int runidjp;
	float jptjp;
	float jneutralfracjp;
	int leadntrkjp;
	int subntrkjp;
	int maxntrkjp;
	int minntrkjp;


	TFile *fmb;
	TTree *treemb;

	int runidmb;
	float jptmb;
	float jneutralfracmb;
	int leadntrkmb;
	int subntrkmb;
	int maxntrkmb;
	int minntrkmb;

	float NeutralFracMax;
	float NeutralFracMin;

	float JetPtCorrMax;

	float mbptcut;


	public:
	ScaleByY(TString yname);
	ScaleByY();
	~ScaleByY();

	bool FillYRatios();	
	double GetYScale(double pt, double nf);
	float GetYScale(float pt, float nf);
	bool Init4Read(TString filenamejp="/home/hep/caines/ly247/Scratch/pp200Y12_jetunderlying/NoTofMatch_FullJet_TransCharged_MatchTrig_ppJP_160811P12id_R06_HadrCorr_VPDcut_161209.root", TString filenamemb="/home/hep/caines/ly247/Scratch/pp200Y12_jetunderlying/NoTofMatch_FullJet_TransCharged_ppMB_160811P12id_R06_HadrCorr_VPDcut_161209.root", TString treename="ResultTree");
	bool WriteY(TString outfilename, bool opt=false);
	

};

#endif



