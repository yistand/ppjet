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
#include "TH2D.h"
#include "TFile.h"
#include "TString.h"
#include "TMatrixD.h"
#include "TLegend.h"
#include "TTree.h"
#include "TFile.h"
#include "TF1.h"

#endif


class TH1D;
class RooUnfoldResponse;
class RooUnfold;
class RooUnfoldErrors;
class RooUnfoldParms;


class Unfold2D {
public:
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

  TString TrigName;

  static const int WNbins = 15;
  static double Wptbins[WNbins+1];

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
  TH2D *hTrain, *hTrainTrue, *hTrainFake, *hTrue, *hMeas, *hFake;
  TH1 *hReco;
  TH1D *hTrainX, *hTrainTrueX, *hTrainFakeX, *hTrueX, *hMeasX, *hRecoX, *hFakeX;
  TH1D *hTrainY, *hTrainTrueY, *hTrainFakeY, *hTrueY, *hMeasY, *hRecoY, *hFakeY;
  TH2D *hCorr;
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

  // Constructors
  Unfold2D (TString name= "EmbedTree.root");
  Unfold2D (const char* name, int argc, const char* const* argv);
  virtual ~Unfold2D() ;

  // Methods and functions
  void  Help();
  void  SetParms (const char* const* argv);
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
  Int_t WriteUnfoldResult();
  Int_t TrainAndTest();
  Int_t Fill4Test(int *Nevents);
  Int_t WriteTest();
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
  TH2D *CorrelationHist (const TMatrixD& cov,const char* name, const char* title,Double_t lo, Double_t hi)   ;

  static TH2D* RebinX2DHisto(TH2D *old);
  static TH2D* RebinAs(TH2D *old, TH2D *model);		// rebin old as the same binining as model
  static TH2D* Rebin2DHisto(TH2D *old, int Nx, double *xbins, int Ny, double *ybins);
  static TH2D *InitWideX2DHisto(TString name, TString title, int Ny, double ylo, double yhi);
  void SetXYname(TString xname, TString yname);
  void SetTrigName(TString tname);

};




//#ifndef NOINLINE
//#include "Unfold2D.cxx"
//#endif

#endif
