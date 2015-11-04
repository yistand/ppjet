//================================================================================
//
//	2015.10.06	Li YI
//	read from underlying event result root files
//	plot each histogram from the root files
//
//================================================================================
#include <iostream>

#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TLatex.h"
#include "TString.h"

using namespace std;


void plot(char *dir = "/home/hep/caines/ly247/Scratch/pp200Y12_jetunderlying/",char *filetemplate = "underlyingevent_%s_R06_LeadJetAngle",char *nametag="JP2") {

	char filetag[500];
	sprintf(filetag,filetemplate,nametag);
	TString filename = Form("%s%s.root",dir,filetag);
	cout<<"Reading "<<filename<<endl;
	TFile *file = new TFile(filename);
	if(file->IsZombie()) {
		cout<<"cannot open "<<filename<<endl;
		return;
	}

	TH1D* LeadJetPt = (TH1D*)file->Get("LeadingJetPt");
	TH1D* SubJetPt = (TH1D*)file->Get("SubLeadingJetPt");

	TProfile* LeadJetNtrkvsLeadJetPt = (TProfile*)file->Get("LeadingJetNtrkvsLeadJetPt");
	TProfile* SubJetNtrkvsLeadJetPt = (TProfile*)file->Get("SubJetNtrkvsLeadJetPt");
	TProfile* TranMaxNtrkvsLeadJetPt = (TProfile*)file->Get("TranMaxNtrkvsLeadJetPt");
	TProfile* TranMinNtrkvsLeadJetPt = (TProfile*)file->Get("TranMinNtrkvsLeadJetPt");
	TProfile* TranNtrkvsLeadJetPt = (TProfile*)file->Get("TranNtrkvsLeadJetPt");


	TProfile* LeadJetPtvsLeadJetPt = (TProfile*)file->Get("LeadingJetPtvsLeadJetPt");
	TProfile* SubJetPtvsLeadJetPt = (TProfile*)file->Get("SubJetPtvsLeadJetPt");
	TProfile* TranMaxPtvsLeadJetPt = (TProfile*)file->Get("TranMaxPtvsLeadJetPt");
	TProfile* TranMinPtvsLeadJetPt = (TProfile*)file->Get("TranMinPtvsLeadJetPt");
	TProfile* TranPtvsLeadJetPt = (TProfile*)file->Get("TranPtvsLeadJetPt");

	TH2D* Spectrum_LeadJetPtvsLeadJetPt = (TH2D*)file->Get("Spectrum_LeadingJetPtvsLeadJetPt");
	TH2D* Spectrum_SubJetPtvsLeadJetPt = (TH2D*)file->Get("Spectrum_SubJetPtvsLeadJetPt");
	TH2D* Spectrum_TranMaxPtvsLeadJetPt = (TH2D*)file->Get("Spectrum_TranMaxPtvsLeadJetPt");
	TH2D* Spectrum_TranMinPtvsLeadJetPt = (TH2D*)file->Get("Spectrum_TranMinPtvsLeadJetPt");
	TH2D* Spectrum_TranPtvsLeadJetPt = (TH2D*)file->Get("Spectrum_TranPtvsLeadJetPt");




	LeadJetPt->SetTitle("Leading Jet Pt");
	LeadJetPt->GetXaxis()->SetTitle("Leading Jet p_{T}");

	SubJetPt->SetTitle("SubLeading Jet Pt");
	SubJetPt->GetXaxis()->SetTitle("SubLeading Jet p_{T}");


	LeadJetNtrkvsLeadJetPt->SetTitle("LeadArea Ntrk vs Leading Jet Pt");
	LeadJetNtrkvsLeadJetPt->GetXaxis()->SetTitle("Leading Jet p_{T}");
	LeadJetNtrkvsLeadJetPt->GetYaxis()->SetTitle("Ntrk in Leading Jet Area");

	SubJetNtrkvsLeadJetPt->SetTitle("SubLeadArea Ntrk vs Leading Jet Pt");
	SubJetNtrkvsLeadJetPt->GetXaxis()->SetTitle("Leading Jet p_{T}");
	SubJetNtrkvsLeadJetPt->GetYaxis()->SetTitle("Ntrk in SubLeading Jet Area");

	TranMaxNtrkvsLeadJetPt->SetTitle("Transverse Max Ntrk vs Leading Jet Pt");
	TranMaxNtrkvsLeadJetPt->GetXaxis()->SetTitle("Leading Jet p_{T}");
	TranMaxNtrkvsLeadJetPt->GetYaxis()->SetTitle("Ntrk in Transverse Max Area");

	TranMinNtrkvsLeadJetPt->SetTitle("Transverse Min Ntrk vs Leading Jet Pt");
	TranMinNtrkvsLeadJetPt->GetXaxis()->SetTitle("Leading Jet p_{T}");
	TranMinNtrkvsLeadJetPt->GetYaxis()->SetTitle("Ntrk in Transverse Min Area");

	TranNtrkvsLeadJetPt->SetTitle("Transverse Ntrk vs Leading Jet Pt");
	TranNtrkvsLeadJetPt->GetXaxis()->SetTitle("Leading Jet p_{T}");
	TranNtrkvsLeadJetPt->GetYaxis()->SetTitle("Ntrk in Transverse  Area");


	LeadJetPtvsLeadJetPt->SetTitle("LeadArea <Pt> vs Leading Jet Pt");
	LeadJetPtvsLeadJetPt->GetXaxis()->SetTitle("Leading Jet p_{T}");
	LeadJetPtvsLeadJetPt->GetYaxis()->SetTitle("<Pt> in Leading Jet Area ");

	SubJetPtvsLeadJetPt->SetTitle("SubLeadArea <Pt> vs Leading Jet Pt");
	SubJetPtvsLeadJetPt->GetXaxis()->SetTitle("Leading Jet p_{T}");
	SubJetPtvsLeadJetPt->GetYaxis()->SetTitle("<Pt> in SubLeading Jet Area ");

	TranMaxPtvsLeadJetPt->SetTitle("Transverse Max <Pt> vs Leading Jet Pt");
	TranMaxPtvsLeadJetPt->GetXaxis()->SetTitle("Leading Jet p_{T}");
	TranMaxPtvsLeadJetPt->GetYaxis()->SetTitle("Transverse Max <Pt> ");

	TranMinPtvsLeadJetPt->SetTitle("Transverse Min <Pt> vs Leading Jet Pt");
	TranMinPtvsLeadJetPt->GetXaxis()->SetTitle("Leading Jet p_{T}");
	TranMinPtvsLeadJetPt->GetYaxis()->SetTitle("Transverse Min <Pt> ");

	TranPtvsLeadJetPt->SetTitle("Transverse <Pt> vs Leading Jet Pt");
	TranPtvsLeadJetPt->GetXaxis()->SetTitle("Leading Jet p_{T}");
	TranPtvsLeadJetPt->GetYaxis()->SetTitle("Transverse  <Pt> ");

	Spectrum_LeadJetPtvsLeadJetPt->SetTitle("LeadArea Pt Spectrum vs Leading Jet Pt");
	Spectrum_LeadJetPtvsLeadJetPt->GetXaxis()->SetTitle("Leading Jet p_{T}");
	Spectrum_LeadJetPtvsLeadJetPt->GetYaxis()->SetTitle("Pt Spectrum in Leading Jet Area");

	Spectrum_SubJetPtvsLeadJetPt->SetTitle("SubLeadArea Pt Spectrum vs Leading Jet Pt");
	Spectrum_SubJetPtvsLeadJetPt->GetXaxis()->SetTitle("Leading Jet p_{T}");
	Spectrum_SubJetPtvsLeadJetPt->GetYaxis()->SetTitle("Pt Spectrum in SubLeading Jet Area");

	Spectrum_TranMaxPtvsLeadJetPt->SetTitle("Transverse Max Pt Spectrum vs Leading Jet Pt");
	Spectrum_TranMaxPtvsLeadJetPt->GetXaxis()->SetTitle("Leading Jet p_{T}");
	Spectrum_TranMaxPtvsLeadJetPt->GetYaxis()->SetTitle("Pt Spectrum in Transverse Max Area");

	Spectrum_TranMinPtvsLeadJetPt->SetTitle("Transverse Min Pt Spectrum vs Leading Jet Pt");
	Spectrum_TranMinPtvsLeadJetPt->GetXaxis()->SetTitle("Leading Jet p_{T}");
	Spectrum_TranMinPtvsLeadJetPt->GetYaxis()->SetTitle("Pt Spectrum in Transverse Min Area");

	Spectrum_TranPtvsLeadJetPt->SetTitle("Transverse Pt Spectrum vs Leading Jet Pt");
	Spectrum_TranPtvsLeadJetPt->GetXaxis()->SetTitle("Leading Jet p_{T}");
	Spectrum_TranPtvsLeadJetPt->GetYaxis()->SetTitle("Pt Spectrum in Transverse  Area");


	TCanvas *c = new TCanvas();
	TLatex *lat = new TLatex(0.09,0.94, nametag);
	//TLatex *lat = new TLatex(0.02,0.94, "JP0+1+2");
	lat->SetNDC();
	lat->SetTextColor(2);

	c->SetLogy();
	LeadJetPt->Draw();
	lat->Draw("same");
	c->SaveAs(Form("LeadingJetPt_%s.png",filetag));
	c->Clear();
	SubJetPt->Draw();
	lat->Draw("same");
	c->SaveAs(Form("SubLeadingJetPt_%s.png",filetag));
	c->Clear();
	c->SetLogy(0);

	LeadJetNtrkvsLeadJetPt->Draw();
	lat->Draw("same");
	c->SaveAs(Form("LeadingJetNtrkvsLeadJetPt_%s.png",filetag));
	c->Clear();
	SubJetNtrkvsLeadJetPt->Draw();
	lat->Draw("same");
	c->SaveAs(Form("SubJetNtrkvsLeadJetPt_%s.png",filetag));
	c->Clear();
	TranMaxNtrkvsLeadJetPt->Draw();
	lat->Draw("same");
	c->SaveAs(Form("TranMaxNtrkvsLeadJetPt_%s.png",filetag));
	c->Clear();
	TranMinNtrkvsLeadJetPt->Draw();
	lat->Draw("same");
	c->SaveAs(Form("TranMinNtrkvsLeadJetPt_%s.png",filetag));
	c->Clear();
	TranNtrkvsLeadJetPt->Draw();
	lat->Draw("same");
	c->SaveAs(Form("TranNtrkvsLeadJetPt_%s.png",filetag));
	c->Clear();


	LeadJetPtvsLeadJetPt->Draw();
	lat->Draw("same");
	c->SaveAs(Form("LeadingJetPtvsLeadJetPt_%s.png",filetag));
	c->Clear();
	SubJetPtvsLeadJetPt->Draw();
	lat->Draw("same");
	c->SaveAs(Form("SubJetPtvsLeadJetPt_%s.png",filetag));
	c->Clear();
	TranMaxPtvsLeadJetPt->Draw();
	lat->Draw("same");
	c->SaveAs(Form("TranMaxPtvsLeadJetPt_%s.png",filetag));
	c->Clear();
	TranMinPtvsLeadJetPt->Draw();
	lat->Draw("same");
	c->SaveAs(Form("TranMinPtvsLeadJetPt_%s.png",filetag));
	c->Clear();
	TranPtvsLeadJetPt->Draw();
	lat->Draw("same");
	c->SaveAs(Form("TranPtvsLeadJetPt_%s.png",filetag));
	c->Clear();

	c->SetLogz();
	Spectrum_LeadJetPtvsLeadJetPt->Draw("col");
	lat->Draw("same");
	c->SaveAs(Form("Spectrum_LeadingJetPtvsLeadJetPt_%s.png",filetag));
	c->Clear();
	Spectrum_SubJetPtvsLeadJetPt->Draw("col");
	lat->Draw("same");
	c->SaveAs(Form("Spectrum_SubJetPtvsLeadJetPt_%s.png",filetag));
	c->Clear();
	Spectrum_TranMaxPtvsLeadJetPt->Draw("col");
	lat->Draw("same");
	c->SaveAs(Form("Spectrum_TranMaxPtvsLeadJetPt_%s.png",filetag));
	c->Clear();
	Spectrum_TranMinPtvsLeadJetPt->Draw("col");
	lat->Draw("same");
	c->SaveAs(Form("Spectrum_TranMinPtvsLeadJetPt_%s.png",filetag));
	c->Clear();
	Spectrum_TranPtvsLeadJetPt->Draw("col");
	lat->Draw("same");
	c->SaveAs(Form("Spectrum_TranPtvsLeadJetPt_%s.png",filetag));
	c->Clear();
	c->SetLogz(0);



}



