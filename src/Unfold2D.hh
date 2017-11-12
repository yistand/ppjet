//=====================================================================-*-C++-*-
//
//	2016.08.23		Li YI
//	try to unfold underlying event activity vs leading jet pt
//
//==============================================================================
//
// File and Version Information:
//      $Id: RooUnfoldTestHarness2D.h 322 2011-10-27 00:23:35Z T.J.Adye $
//
// Description:
//      Test Harness class for the RooUnfold package using 2D toy MC.
//      Inherits from RooUnfoldTestHarness.
//
// Author: Tim Adye <T.J.Adye@rl.ac.uk>
//
//==============================================================================

#ifndef UNFOLD2D_HH
#define UNFOLD2D_HH

#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TH1.h"
#include "TH2.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TFile.h"
#include "TString.h"
#include "TMatrixD.h"
#include "TLegend.h"
#include "TTree.h"
#include "TFile.h"
#include "TF1.h"

#include <vector>

#endif


class TH1D;
class RooUnfoldResponse;
class RooUnfold;
class RooUnfoldErrors;
class RooUnfoldParms;


class Unfold2D {
private:
  // Parameters
  Int_t    method;
  Int_t    ntx, nmx;
  Int_t    nty, nmy;
  Double_t xlo, xhi;
  Double_t ylo, yhi;
  Int_t    overflow;
  Int_t    verbose;		
  Int_t    doerror;
  Int_t    regparm;
  bool   WIDEBIN;
  bool	 flagjetweight;
  bool	flagscaley;

  // MC root file options
  bool OPT_TPCSYS;
  bool OPT_TPCSYSPM;		// true for plus, false for minus 
  bool OPT_TPCSYSABS;		// true use absolute, false use relative 5%

  bool OPT_BEMCSYS;
  bool OPT_BEMCSYSPM;		// true for plus, false for minus  4%

  // ==1: RC>0.2, Mc>0.5, in order to compare with other experiments. ==2: RC>0.5, Mc>0.5
  int OPT_RC02MC05;

  // MIP instead of default 100% hadronic correction
  bool OPT_MIP;

  // weight Rc Vz to be same as real data
  bool OPT_DORCVZWEIGHT;

  // Unfolding Options
  bool NoFakeOpt;		// ==1: no fake, ==0: fake included (default)
  bool NoLossOpt;		// ==1: no eff, ==0: eff included (default)

  TString TrigName;

  TString TranCharge;

  Int_t ExcludeOpt;

  Int_t CHANGEPRIOR;

  //static const int WNbins = 13;
  static const int WNbins = 12;
  static double Wptbins[WNbins+1];

  static const int MAXARRAY = 500;

  int FlagDefaultCalled;

  RooUnfoldResponse* response;
  RooUnfold*         unfold;


  // Training Data
  // --------------
  TString inputname;
  TFile *ftrain;
  // X-axis and y-axis names for histogram
  TString XvariableName;
  TString YvariableName;
  // Histogram
  TH2F *hTrain, *hTrainTrue, *hTrainFake, *hTrue, *hMeas, *hFake;
  TH1 *hReco;
  TH1D *hTrainX, *hTrainTrueX, *hTrainFakeX, *hTrueX, *hMeasX, *hRecoX, *hFakeX;
  TH1D *hTrainY, *hTrainTrueY, *hTrainFakeY, *hTrueY, *hMeasY, *hRecoY, *hFakeY;
  TH2F *hCorr;		
  TH2F *hXRcVsMc, *hYRcVsMc;
  TProfile *pfxTrain; 
  TProfile *pfxTrainTrue; 
  TProfile *pfxMeas; 
  TH1F *hpriorratio;	// if need to change prior, use ratio of measured JP data/pythia MB MC leading jet pt distribution

  TH1F *hMeasVz;	// save measured Vz distribution. If NF weight is applied, it will be the one after NF weight. 
  TH1F *hTrainVz;	// save Train Rc Vz distribution. If Vz weight is applied, it will be the one after Vz weight. 

  // Tree
  TTree *tree;
  // Tree branch
  int runid;
  int eventid;
  float McJet, McPart;
  float RcJet, RcPart;
  int flag;		// 1: true event with Mc and Rc matched, 0: true Mc, no Rc, -1: no Mc, fake Rc
  double weight;	// weighted by cross section / number of the events in that pT bin.

  

  // Output file for training response matrix
  TFile *ftout;

public:
  // Constructors
  Unfold2D (TString name= "EmbedTree.root");
  Unfold2D (const char* name, int argc, const char* const* argv);
  virtual ~Unfold2D() ;

  // Methods and functions
  void  Help();
  void  SetParms (const char* const* argv);
  void  SetNoFake(bool val=true) {NoFakeOpt = val;}
  void  SetNoLoss(bool val=true) {NoLossOpt = val;}
  bool  GetNoFakeStatus() {return NoFakeOpt;}
  bool  GetNoEffStatus() {return NoLossOpt;}
  int   Float2Int(float aFloat);
  void  SetParms (float* args);
  void  SetDefaultParms ();
  void  PrintParms ();
  Int_t FillbyXsec4Train ();
  Int_t FillbyXsec4Train (int *Nevents);
  Int_t Fill4Train ();
  Int_t ReadResponseMatrix();
  Int_t Train();
  void TrainResults();
  Int_t WriteTrain();
  Int_t Fill4Unfold ();
  Int_t ReadMeasHist4Unfold();
  void  Results();
  Int_t Unfold();
  Int_t WriteHist4Unfold();	// for Debug
  Int_t WriteUnfoldResult();
  Int_t TrainAndTest();
  Int_t Fill4Test(int *Nevents);
  Int_t WriteTest();
  float Weight2ChargePrior(float jetpt);
  float Weight2ChargePrior(TF1 *f1, float jetpt);
  virtual void  Reset();
  virtual void  Init();
  virtual Int_t CheckParms();
  static TH1D* ProjectionX (const TH1* h, const char* name="_px", const char* title=0, Option_t* opt="");		// By declaring a function member as static, you make it independent of any particular object of the class. A static member function can be called even if no objects of the class exist and the static functions are accessed using only the class name and the scope resolution operator ::
  static TH1D* ProjectionY (const TH1* h, const char* name="_py", const char* title=0, Option_t* opt="");
  static TProfile* ProfileX (const TH1* h, const char* name="_pfx", const char* title=0, Option_t* opt="");
  static TH1D* TProfile2TH1D(TProfile *pf);
  static void     setmax   (TH1* h, const TH1* h1= 0, const TH1* h2= 0, const TH1* h3= 0,
                             const TH1* h4= 0, const TH1* h5= 0, const TH1* h6= 0);
  static  void     Legend (TLegend*& legend, TH1* truth, TH1* fake, TH1* meas, TH1* reco= 0, TF1* ff=0, TF1* tf=0);
  TH2F *CorrelationHist (const TMatrixD& cov,const char* name, const char* title,Double_t lo, Double_t hi)   ;

  static TH2F* RebinX2DHisto(TH2F *old);
  static TH2F* RebinAs(TH2F *old, TH2F *model);		// rebin old as the same binining as model
  static TH2F* Rebin2DHisto(TH2F *old, int Nx, double *xbins, int Ny, double *ybins);
  static TH2F *InitWideX2DHisto(TString name, TString title, int Ny, double ylo, double yhi);
  static TH2F *InitWideX2DHisto(TString name, TString title, int Ny, double *ybins);
  TH2F *InitWideXY2DHisto(TString name, TString title);
  TH2F *InitWideXX2DHisto(TString name, TString title);
  TH2F *InitWideYY2DHisto(TString name, TString title);
  void SetXYname(TString xname, TString yname);
  void SetTrigName(TString tname);
  void SetTranCharge(TString tname);
  void SetExcludeJPTrig(Int_t exclude);
  void SetChangePrior(Int_t change);
  void SetTpcSys(int opt_tpcsys, int opt_tpcsyspm, int opt_tpcsysabs);	// TPC tracking 5% efficiency study (opt_tpcsys==1), use absolute (opt_tpcsysabs==1) or relative plus(opt_tpcsyspm==1) or minus 5% 
  void SetBemcSys(int opt_bemcsys, int opt_bemcsyspm);	// BEMC gain 4% uncertainty study (opt_bemcsys==1), plus(opt_bemcsyspm==1) or minus 4% 
  void SetRc02Mc05(int opt_rc02mc05);	// ==1 : RC>0.2, Mc>0.5, in order to compare with other experiments only unfolded to pt>0.5.  ==2 : RC>0.5, Mc>0.5
  void SetMIP(int opt_mip);	// read MIP root file instead of 100% hadronic correction
  void SetDoRcVzWeight(int opt_dorcvzweight);	// weight Rc Vz to be same as real data
  int Find(std::vector<int> vec, int val);

  void GetTpcErrString(TString &STpcSys);
  

};




//#ifndef NOINLINE
//#include "Unfold2D.cxx"
//#endif

#endif
