#ifndef ROOT_REWEIGHTBYNF
#define ROOT_REWEIGHTBYNF

#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"

class ReWeightByNF{
	
	TH2F *hw;
	
	TH2F *hjp;
	TH2F *hmb;

	TH1D *htmprx;

	TFile *fjp;
	TTree *treejp;

	int runidjp;
	float jptjp;
	float jneutralfracjp;

	TFile *fmb;
	TTree *treemb;

	int runidmb;
	float jptmb;
	float jneutralfracmb;

	float NeutralFracMax;
	float NeutralFracMin;

	float JetPtCorrMax;

	float mbptcut;


	public:
	ReWeightByNF();
	~ReWeightByNF();

	bool FillNFRatios();	
	double GetNFWeight(double pt, double nf);
	float GetNFWeight(float pt, float nf);
	bool AdjustNevts(TH2D *h2, TH1D *hx);
	bool AdjustNevts(TH2F *h2, TH1D *hx);
	//bool Init4Read(TString filenamejp="/home/hep/caines/ly247/Scratch/pp200Y12_jetunderlying/NoTofMatch_FullJet_TransCharged_MatchTrig_ppJP_160811P12id_R06_HadrCorr_VPDcut_161209.root", TString filenamemb="/home/hep/caines/ly247/Scratch/pp200Y12_jetunderlying/NoTofMatch_FullJet_TransCharged_ppMB_160811P12id_R06_HadrCorr_VPDcut_161209.root", TString treename="ResultTree");
	bool Init4Read(TString filenamejp="/home/hep/caines/ly247/Scratch/pp200Y12_jetunderlying/NoTofMatch_FullJet_TransCharged_MatchTrig_ppJP_160811P12id_R06_HadrCorr_VPDcut_170418.root", TString filenamemb="/home/hep/caines/ly247/Scratch/pp200Y12_jetunderlying/NoTofMatch_FullJet_TransCharged_ppMB_160811P12id_R06_HadrCorr_VPDcut_170418.root", TString treename="ResultTree");
	bool WriteNF(TString outfilename, bool opt=false);
	bool ReadNFfile(TString infilename);
	

};

#endif



