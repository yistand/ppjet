//====================================================================================================
//
//		2016.07.24	Li Yi
//		Merge Unfolding matrix from each pT bin by their cross section
//
//
//		2016.07.25	Li Yi
//		pt2_3 has one entry with extremely high Mc pt. This may cause trouble when merging
//		because low pT bin will have huge weight, therefore we delete this bin
//
//====================================================================================================

//#include "HCCrossSectionPerpT.h"
#include "CrossSectionPerpT.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TColor.h"

#include <iostream>


int main() {

	TString flagLeadSub="LeadOrSub"; // "Lead", "LeadOrSub"

	TH2D *Cov;			// When Rc matched Leading or SubLeading Mc Jet  (x: Rc; Y: Mc);
	TH1D *Mc;			// Mc Leading jet distribution no matter matched or not
	TH1D *Rc;			// Rc Leading jet distribution no matched matched to Mc or not

	TH1D *summc;			// not weighted Mc Leading jet distribution no matter matched or not
	TH1D *sumrc;			// not weighted Rc Leading jet distribution no matched matched to Mc or not

	TH2D *tmpcov;
	TH1D *tmpmc;
	TH1D *tmprc;

	int NpT = NUMBEROFPT;		// defined in include/CrossSectionPerpT.h  

	//int color[11] = {kBlack, kBlue, kBlue-7, kViolet+1, kViolet-1, kMagenta-7, kMagenta, kPink+10, kPink-9, kRed-7, kRed};
	int color[11] = {1,2,kGreen+2,4,kOrange,6,kCyan+1,8,9,12,41};
	TCanvas *cperptMc = new TCanvas("cperptMc");
	cperptMc->SetLogy();
	TCanvas *cbyxsecMc = new TCanvas("cbyxsecMc");
	cbyxsecMc->SetLogy();

	TCanvas *cperptRc = new TCanvas("cperptRc");
	cperptRc->SetLogy();
	TCanvas *cbyxsecRc = new TCanvas("cbyxsecRc");
	cbyxsecRc->SetLogy();

	TFile *fin;
	for(int i = 0; i<NpT ; i++) {
		TString ifilename = TString("~/Scratch/embedPythia/pt")+TString(PTBINS[i])+TString("_JetMcVsEmbedMatchTrig_Online_nobbc.root");
		//TString ifilename = TString("~/Scratch/embedPythia/pt")+TString(PTBINS[i])+TString("_JetMcVsEmbedMatchTrig.root");
		//TString ifilename = TString("~/Scratch/embedPythia/pt")+TString(PTBINS[i])+TString("_JetMcVsEmbed.root");
		//TString ifilename = TString("~/Scratch/embedPythia/HCpt")+TString(PTBINS[i])+TString("_JetMcVsEmbedMatchTrig.root");
		fin = new TFile(ifilename);
		if(fin->IsOpen()&&!fin->IsZombie()) {
			std::cout<<"reading "<<fin->GetName()<<std::endl;
			tmpcov = (TH2D*)fin->Get("CovMcMatched"+flagLeadSub+"JetVsRcJet");
			tmpmc = (TH1D*)fin->Get("McLeadJetPt");
			tmprc = (TH1D*)fin->Get("RcLeadJetPt");

			if(!tmpcov || !tmpmc || !tmprc) std::cout<<"cannot find histogram"<<std::endl;
			//std::cout<<tmpcov->Integral()<<"\t"<<tmpmc->Integral()<<"\t"<<tmprc->Integral()<<std::endl;
			
			if(i==0) {
				Cov = (TH2D*)tmpcov->Clone("CovMatrix");
				Cov->Sumw2();
				Cov->Scale(XSEC[0]/NUMBEROFEVENT[0]);
				Mc = (TH1D*)tmpmc->Clone("Mc");
				Mc->Sumw2();
				Mc->Scale(XSEC[0]/NUMBEROFEVENT[0]);
				Rc = (TH1D*)tmprc->Clone("Rc");
				Rc->Sumw2();
				Rc->Scale(XSEC[0]/NUMBEROFEVENT[0]);
				std::cout<<PTBINS[i]<<" Mc integral="<<tmpmc->Integral()<<"\tRc integral="<<tmprc->Integral()<<"\tratio="<<1.*tmpmc->Integral()/tmprc->Integral()<<std::endl;
				std::cout<<"Mc integral="<<Mc->Integral()<<"\tRc integral="<<Rc->Integral()<<std::endl;


				Cov->GetXaxis()->SetTitle("Detector-level Leading Jet p_{T}");
				Cov->GetYaxis()->SetTitle("Particle-level Leading Jet p_{T}");
				Mc->GetXaxis()->SetTitle("Particle-level Leading Jet p_{T}");
				Rc->GetXaxis()->SetTitle("Detector-level Leading Jet p_{T}");


				summc = (TH1D*)tmpmc->Clone("summc");
				summc->Sumw2();
				sumrc = (TH1D*)tmprc->Clone("sumrc");
				sumrc->Sumw2();
			
				cperptMc->cd();
				tmpmc->SetLineWidth(2);
				tmpmc->SetLineColor(color[i]);
				tmpmc->DrawClone("h");

				cbyxsecMc->cd();
				tmpmc->Scale(XSEC[0]/NUMBEROFEVENT[0]);
				tmpmc->DrawClone("h");
				tmpmc->Scale(1/(XSEC[0]/NUMBEROFEVENT[0]));


				cperptRc->cd();
				tmprc->SetLineWidth(2);
				tmprc->SetLineColor(color[i]);
				tmprc->DrawClone("h");
				//for(int jj = 0; jj<tmprc->GetNbinsX();jj++) std::cout<<tmprc->GetBinCenter(jj+1)<<"\t"<<tmprc->GetBinContent(jj+1)<<std::endl;

				cbyxsecRc->cd();
				tmprc->Scale(XSEC[0]/NUMBEROFEVENT[0]);
				tmprc->DrawClone("h");
				tmprc->Scale(1/(XSEC[0]/NUMBEROFEVENT[0]));
			}
			else {

				Cov->Add(tmpcov,XSEC[i]/NUMBEROFEVENT[i]);
				Mc->Add(tmpmc,XSEC[i]/NUMBEROFEVENT[i]);
				Rc->Add(tmprc,XSEC[i]/NUMBEROFEVENT[i]);

				std::cout<<PTBINS[i]<<" Mc integral="<<tmpmc->Integral()<<"\tRc integral="<<tmprc->Integral()<<"\tratio="<<1.*tmpmc->Integral()/tmprc->Integral()<<std::endl;
				std::cout<<"Mc integral="<<Mc->Integral()<<"\tRc integral="<<Rc->Integral()<<std::endl;

				summc->Add(tmpmc);
				sumrc->Add(tmprc);

				cperptMc->cd();
				tmpmc->SetLineWidth(2);
				tmpmc->SetLineColor(color[i]);
				tmpmc->DrawClone("hsame");

				cbyxsecMc->cd();
				tmpmc->Scale(XSEC[i]/NUMBEROFEVENT[i]);
				tmpmc->DrawClone("hsame");
				tmpmc->Scale(1/(XSEC[i]/NUMBEROFEVENT[i]));


				cperptRc->cd();
				tmprc->SetLineWidth(2);
				tmprc->SetLineColor(color[i]);
				tmprc->DrawClone("hsame");

				//for(int jj = 0; jj<tmprc->GetNbinsX();jj++) std::cout<<tmprc->GetBinCenter(jj+1)<<"\t"<<tmprc->GetBinContent(jj+1)<<std::endl;

				cbyxsecRc->cd();
				tmprc->Scale(XSEC[i]/NUMBEROFEVENT[i]);
				tmprc->DrawClone("hsame");
				tmprc->Scale(1/(XSEC[i]/NUMBEROFEVENT[i]));
			}	
			//fin->Close();		// histogram is lost if closed

		}
	}
	//fin->Delete();		// if delete, the last histogram is lost..


	TFile *fout = new TFile("UnfoldMatrx"+flagLeadSub+"_Online_nobbc.root","RECREATE");
	//TFile *fout = new TFile("HCUnfoldMatrx"+flagLeadSub+".root","RECREATE");
	std::cout<<"Write to "<<fout->GetName()<<std::endl;
	fout->cd();
	Cov->Write();
	Mc->Write();
	Rc->Write();

	summc->Write();
	sumrc->Write();
	
	cperptMc->Write();
	cbyxsecMc->Write();
	
	cperptRc->Write();
	cbyxsecRc->Write();
	
	fout->Close();

	return 1;
}

