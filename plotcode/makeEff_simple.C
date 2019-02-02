#include "TCanvas.h"
#include "TSystem.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TString.h"
#include "TMath.h"
#include "math.h"
#include <iostream>

using namespace std;

int ZdcArray(float zdc) {

  int array = 0;
  array = int(zdc/1000)-3;

  if(array>=10) array = 9;
  if(array<0) array = 0;

  return array; 

}


int GlobalArray(float globaltrk) {

  int array = 0;
  array = int(globaltrk/100);

  if(array>=10) array = 9;

  return array; 

}


void makeEff_simple(const char* particle="Piplus", int energy=200, TString dependOn="RefMult"){	
//	dependOn: RefMult, GlobalTr, Zdc
  if( (!dependOn.EqualTo("RefMult",TString::kIgnoreCase)) && (!dependOn.EqualTo("GlobalTr",TString::kIgnoreCase)) && (!dependOn.EqualTo("Zdc",TString::kIgnoreCase)) ) {
    cout<<"ERR!! call makeEff_simple(const char* particle, int energy, TString dependOn): dependOn should be \"RefMult\", \"GlobalTr\" or \"Zdc\"."<<endl;
    return;
  }

  static const Double_t pi = TMath::Pi();

  gStyle->SetOptStat(111111);
  gStyle->SetPalette(1);
  //gStyle->SetFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetFrameBorderMode(0);
  //gStyle->SetPadTopMargin(0.05);
  //gStyle->SetPadRightMargin(0.05);
  //gStyle->SetPadBottomMargin(0.05);
  gStyle->SetPadLeftMargin(0.11);
  TF1* feff = new TF1("feff","[0]*pow(1.-[1]*(1.-[2])*x*x,1./(1.-[2]))",0.10,10.0);
  feff->SetParameters(63.2718,4.55971,1.1784);
  
  int Nbins = 100;
  float max_z = 30;
  float min_pt = 0.;
  float max_pt = 20;
  float min_fitpts = 25;
  float min_fitpts_nposs = .52;
  float max_eta = 1.0;
  float max_dca = 1.0;
  float PID = 0;
  char cbuff[10];
  char buffer[100];

  sprintf(cbuff, "%d",energy);
  TString NRG = TString(cbuff);
  
  cout << "particle = " << particle << endl;
  cout << "energy = " << NRG << endl;
  
  TString out_file, in_file;
  TString prefix; 
  TString postfix; 
  prefix = TString(".");
  postfix = TString(".");
  double mass;
  if(TString(particle) == TString("Piplus")){
    in_file = prefix + TString("/SinglePiplusNT_Embed_") + NRG + TString("GeV.root");
    out_file = postfix + TString("/") + dependOn + TString("DepPiplus") + NRG + TString("GeV.root");
    mass = .13975;
    PID = 8;
  }
  else if(TString(particle) == TString("Piminus")){
    in_file = prefix + TString("/SinglePiminusNT_Embed_") + NRG + TString("GeV.root");
    out_file = postfix + TString("/") + dependOn + TString("DepPiminus") + NRG + TString("GeV.root");
    mass = .13975;
    PID = 9;
  }
  else if(TString(particle) == TString("Kplus")){
    in_file = prefix + TString("/SingleKplusNT_Embed_") + NRG + TString("GeV.root");
    out_file = postfix + TString("/") + dependOn + TString("DepKplus") + NRG + TString("GeV.root");
    mass = .493677;
    PID = 11;
  }
  else if(TString(particle) == TString("Kminus")){
    in_file = prefix + TString("/SingleKminusNT_Embed_") + NRG + TString("GeV.root");
    out_file = postfix + TString("/") + dependOn + TString("DepKminus") + NRG + TString("GeV.root");
    mass = .493677;
    PID = 12;
  }
  else if(TString(particle) == TString("Proton")){
    in_file = prefix + TString("/SingleProtonNT_Embed_") + NRG + TString("GeV.root");
    out_file = postfix + TString("/") + dependOn + TString("DepProton") + NRG + TString("GeV.root");
    mass = .93827;
    PID = 14;
  }
  else if(TString(particle) == TString("AntiProton")){
    in_file = prefix + TString("/SingleAntiProtonNT_Embed_") + NRG + TString("GeV.root");
    out_file = postfix + TString("/") + dependOn + TString("DepAntiProton") + NRG + TString("GeV.root");
    mass = .93827;
    PID = 15;
  }
  

//------------------------- Input ---------------------------
  cout << "opening input file: " << in_file << endl;
  TFile *f = new TFile(in_file, "READ");
  TNtuple *MatchedPairs = (TNtuple*) f->Get("MatchedPairs_NT");
  TNtuple *McTrack = (TNtuple*) f->Get("McTrack_NT");
  
  float Dedx, RefMult, RefMultCorrected, CentralityWeight, Centrality16, VertexX, VertexY, VertexZ, PtMc, PzMc, EtaMc, PhiMc, PtPr, PzPr, EtaPr, PhiPr, DcaGl, DcaZGl, DcaXYGl, Flag, FitPts, DedxPts, AllPts, NPossible, ParentGeantId, GeantId,mErrP, GlobalTr, ZdcRate, BbcRate;
  MatchedPairs->SetBranchAddress("Dedx", &Dedx);
  MatchedPairs->SetBranchAddress("RefMult", &RefMult);
  MatchedPairs->SetBranchAddress("RefMultCorrected", &RefMultCorrected);
  MatchedPairs->SetBranchAddress("CentralityWeight", &CentralityWeight);
  MatchedPairs->SetBranchAddress("Centrality16", &Centrality16);
  MatchedPairs->SetBranchAddress("VertexX", &VertexX);
  MatchedPairs->SetBranchAddress("VertexY", &VertexY);
  MatchedPairs->SetBranchAddress("VertexZ", &VertexZ);
  MatchedPairs->SetBranchAddress("PtMc", &PtMc);
  MatchedPairs->SetBranchAddress("PzMc", &PzMc);
  MatchedPairs->SetBranchAddress("EtaMc", &EtaMc);
  MatchedPairs->SetBranchAddress("PhiMc", &PhiMc);
  MatchedPairs->SetBranchAddress("PtPr", &PtPr);
  MatchedPairs->SetBranchAddress("mErrP", &mErrP);
  //MatchedPairs->SetBranchAddress("PzPr", &PzPr);
  MatchedPairs->SetBranchAddress("EtaPr", &EtaPr);
  MatchedPairs->SetBranchAddress("PhiPr", &PhiPr);
  MatchedPairs->SetBranchAddress("DcaGl", &DcaGl);
  MatchedPairs->SetBranchAddress("DcaZGl", &DcaZGl);
  MatchedPairs->SetBranchAddress("DcaXYGl", &DcaXYGl);
  MatchedPairs->SetBranchAddress("Flag", &Flag);
  MatchedPairs->SetBranchAddress("FitPts", &FitPts);
  MatchedPairs->SetBranchAddress("DedxPts", &DedxPts);
  MatchedPairs->SetBranchAddress("AllPts", &AllPts);
  MatchedPairs->SetBranchAddress("NPossible", &NPossible);
  MatchedPairs->SetBranchAddress("ParentGeantId", &ParentGeantId);
  MatchedPairs->SetBranchAddress("GeantId", &GeantId);
  MatchedPairs->SetBranchAddress("GlobalTr", &GlobalTr);
  MatchedPairs->SetBranchAddress("ZdcRate", &ZdcRate);
  MatchedPairs->SetBranchAddress("BbcRate", &BbcRate);
  
  float pRefMult, pRefMultCorrected, pCentralityWeight, pCentrality16, pVertexX, pVertexY, pVertexZ, pPtMc, pPzMc, pEtaMc, pPhiMc, pParentGeantId, pGeantId, pGlobalTr, pZdcRate, pBbcRate;
  McTrack->SetBranchAddress("RefMult", &pRefMult);
  McTrack->SetBranchAddress("RefMultCorrected", &pRefMultCorrected);
  McTrack->SetBranchAddress("CentralityWeight", &pCentralityWeight);
  McTrack->SetBranchAddress("Centrality16", &pCentrality16);
  McTrack->SetBranchAddress("VertexX", &pVertexX);
  McTrack->SetBranchAddress("VertexY", &pVertexY);
  McTrack->SetBranchAddress("VertexZ", &pVertexZ);
  McTrack->SetBranchAddress("PtMc", &pPtMc);
  McTrack->SetBranchAddress("PzMc", &pPzMc);
  McTrack->SetBranchAddress("EtaMc", &pEtaMc);
  McTrack->SetBranchAddress("PhiMc", &pPhiMc);
  McTrack->SetBranchAddress("ParentGeantId", &pParentGeantId);
  McTrack->SetBranchAddress("GeantId", &pGeantId);
  McTrack->SetBranchAddress("GlobalTr", &pGlobalTr);
  McTrack->SetBranchAddress("ZdcRate", &pZdcRate);
  McTrack->SetBranchAddress("BbcRate", &pBbcRate);
  
//------------------------- Output ---------------------------
  cout << "opening output file: " << out_file << endl;
  TFile *f_out = new TFile(out_file, "RECREATE");
  f_out->cd();
  
  const int Nmult = 10;			// multiplicity range 

  TH1D* hpt_mc[Nmult];
  TH1D* hpt_re[Nmult];
  TH1D* hpt[Nmult];

  TH1D *hallpt_mc;
  TH1D *hallpt_rc;
  TH1D *hallpt;

  for(int i = 0; i<Nmult; i++) {
	hpt_mc[i] = new TH1D(Form("hpt_mc%d",i), Form("MC Tracks vs pt for %s == %d", dependOn.Data(), i), Nbins, 0., 10.);
	hpt_mc[i]->Sumw2();
	hpt_re[i] = new TH1D(Form("hpt_re%d",i), Form("Matched Tracks vs pt for %s == %d", dependOn.Data(), i), Nbins, 0., 10.);
	hpt_re[i]->Sumw2();
  }
  hallpt_mc = new TH1D(Form("hallpt_mc"), Form("MC Tracks vs pt"), Nbins, 0., 10.);
  hallpt_re = new TH1D(Form("hallpt_re"), Form("Matched Tracks vs pt"), Nbins, 0., 10.);

  cout << "histograms created" << endl;
  
//------------------------- Loop ---------------------------
  for(int i = 0; i<McTrack->GetEntries(); i++){
    McTrack->GetEntry(i);
 
    if( pParentGeantId!=0 ) continue;
   
    if( pGeantId!=PID ) continue;

    if( fabs(pVertexZ) > max_z ) continue;

    if( fabs(pEtaMc) > max_eta ) continue;
    if( fabs(pPtMc) < min_pt || fabs(pPtMc) > max_pt) continue;

      if(TString(dependOn) == TString("RefMult")) { 
        if(pRefMult<Nmult) { int j = floor(pRefMult); hpt_mc[j]->Fill(pPtMc); }
      }
      if(TString(dependOn) == TString("GlobalTr")) { 
        //if(pGlobalTr<Nmult) { int j = floor(pGlobalTr); hpt_mc[j]->Fill(pPtMc); }
        int j = GlobalArray(pGlobalTr);
        hpt_mc[j]->Fill(pPtMc);
      }
      if(TString(dependOn) == TString("Zdc")) { 
        int j = ZdcArray(pZdcRate);
	//cout<<j<<endl;
        hpt_mc[j]->Fill(pPtMc);
      }

	hallpt_mc->Fill(pPtMc);
  }
  

  for(int i = 0; i <MatchedPairs->GetEntries(); i++){
    MatchedPairs->GetEntry(i);
    if( ParentGeantId!=0 ) continue;
    if( GeantId!=PID ) continue;
    if( Flag < 0 ) continue;
    if( fabs(VertexZ) > max_z ) continue;
    if( fabs(PtMc)< min_pt || fabs(PtMc)>max_pt) continue;
    if( FitPts < min_fitpts ) continue;
    if( FitPts/NPossible < min_fitpts_nposs ) continue;
    if(fabs(DcaGl) > max_dca ) continue;//change DcaGl-->GcaXYGl
    if( fabs(EtaMc) > max_eta ) continue;
    
    if(TString(dependOn) == TString("RefMult")) { 
      if(RefMult<Nmult) { int j = floor(RefMult); hpt_re[j]->Fill(PtMc); }
    }
    if(TString(dependOn) == TString("GlobalTr")) { 
      //if(GlobalTr<Nmult) { int j = floor(GlobalTr); hpt_re[j]->Fill(PtMc); }
      int j = GlobalArray(GlobalTr);
      hpt_re[j]->Fill(PtMc);
    }
    if(TString(dependOn) == TString("Zdc")) { 
      int j = ZdcArray(ZdcRate);
      //cout<<j<<endl;
      hpt_re[j]->Fill(PtMc);
    }
    
    hallpt_re->Fill(PtMc);
  }

  for(int i = 0; i<Nmult; i++) {
    hpt[i] = (TH1D*)hpt_re[i]->Clone(Form("hpt%d",i));
    hpt[i]->SetTitle(Form("eff vs pt for %s no Cuts == %d",dependOn.Data(),i));
    hpt[i]->Divide(hpt_mc[i]);
  }
  hallpt = (TH1D*)hallpt_re->Clone(Form("hrptall"));
  hallpt->SetTitle(Form("eff vs pt"));
  hallpt->Divide(hallpt_mc);

  f_out->cd();
  f_out->Write();
  for(int i = 0; i<Nmult; i++) {
    hpt_mc[i]->Write();
    hpt_re[i]->Write();
    hpt[i]->Write();
  }
  hallpt_mc->Write();
  hallpt_re->Write();
  hallpt->Write();
  f_out->Close();

  f->Close();

}
