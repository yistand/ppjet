#include "src/UnfoldJetpT2.cxx"


void backfolding(){


	TFile *fin = new TFile("testBayes4_W.root");
	TH1D *unfold = (TH1D*)fin->Get("unfold");
	TH1D *data = (TH1D*)fin->Get("WLeadJetPt_data");
	TH1D *mc = (TH1D*)fin->Get("Wtruth");
	TH1D *rc = (TH1D*)fin->Get("Wmeasured");
	TH2D *cov = (TH2D*)fin->Get("Wcov");
	
	//TH1D *unfold = (TH1D*)fin->Get("unfold");
	//TH1D *data = (TH1D*)fin->Get("LeadJetPt_data");
	//TH1D *mc = (TH1D*)fin->Get("Mc");
	//TH1D *rc = (TH1D*)fin->Get("Rc");
	//TH2D *cov = (TH2D*)fin->Get("CovMatrix");

	//TH1D *unfold = (TH1D*)fin->Get("unfold_train");
	//TH1D *data = (TH1D*)fin->Get("Rc");
	//TH1D *mc = (TH1D*)fin->Get("Mc");
	//TH1D *rc = (TH1D*)fin->Get("Rc");
	//TH2D *cov = (TH2D*)fin->Get("CovMatrix");
	cov->Sumw2();


	TH1D *prorc = (TH1D*)cov->ProjectionX("prorc");
	TH1D *promc = (TH1D*)cov->ProjectionY("promc");
	promc->Sumw2();
	prorc->Sumw2();


	cout<<"get eff and fake"<<endl;


	TH1D *eff = (TH1D*)promc->Clone("eff");
	eff->Divide(mc);
	eff->GetXaxis()->SetTitle("Particle-level Leading Jet p_{T}");
	eff->GetYaxis()->SetTitle("Matched Particle-level Efficiency");


	TH1D *fake = (TH1D*)rc->Clone("fakedist");
	fake->Sumw2();
	fake->Add(prorc,-1);


	fake->GetXaxis()->SetTitle("Detector-level Leading Jet p_{T}");
	fake->GetYaxis()->SetTitle("Unmatched Detector-level Leading Jets");

	TH1D *effunfold = (TH1D*)unfold->Clone("effunfold");
	effunfold->Multiply(eff);

	TH1D *backfold = (TH1D*)data->Clone("backfold");
	backfold->Reset();
	backfold->GetXaxis()->SetTitle("Detector-level Leading Jet p_{T}");
	backfold->GetYaxis()->SetTitle("N_{Leading jets}/dp_{T}");
	TH1D *backfold_nofake = (TH1D*)backfold->Clone("backfold_nofake");
	
	TH1D *hprof;
	for(int i = 1; i<=effunfold->GetNbinsX(); i++) {
		hprob=(TH1D*) cov-> ProjectionX(Form("hprob%d",i),i,i);
		for(int j = 1; j<=effunfold->GetBinContent(i); j++) {
			double r = hprob->GetRandom();
			backfold_nofake->Fill(r);
			backfold->Fill(r);
		}
	}



	//float nfake = 2.57883e+06; // GetnFake(mc, cov, rc, effunfold) ;		fine binning
	float nfake = 2.11342e+06; // GetnFake(mc, cov, rc, effunfold) ;
	cout<<"nfake = "<<nfake<<endl;
	for(int j = 1; j<=nfake; j++) {
		double r = fake->GetRandom();
		backfold->Fill(r);
	}
	float norm = nfake/fake->Integral();
	fake->Scale(norm);

	
	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	TCanvas *c1 = new TCanvas();
	eff->Draw();

	TCanvas *c2 = new TCanvas();
	c2->SetLogy();
	fake->Draw("hist");

	TCanvas *c2r = new TCanvas();
	TH1D *rfake = (TH1D*)fake->Clone("rfake");
	rfake->Divide(data);
	rfake->GetYaxis()->SetTitle("Umatched Detector-level fraction");
	rfake->SetMaximum(1);
	rfake->SetMinimum(0);
	rfake->Draw();

	TCanvas *c3 = new TCanvas();
	c3->SetLogy();
	//backfold->Scale(1,"width");
	//data->Scale(1,"width");
	backfold->SetLineColor(2);
	backfold->SetMarkerColor(2);
	backfold->SetMarkerStyle(8);
	backfold->Draw();
	data->SetMarkerStyle(8);
	data->Draw("same");
	TLegend *leg = new TLegend(0.6,0.7,0.85,0.85);
	leg->SetFillColor(0);
	leg->SetTextFont(42);
	leg->AddEntry(data,"data");
	leg->AddEntry(backfold,"backfold");
	leg->Draw();
	
	TCanvas *c4 = new TCanvas();
	TH1D *ratio = backfold->Clone("ratio");
	ratio->Divide(data);
	ratio->SetMarkerStyle(8);
	ratio->SetMarkerColor(ratio->GetLineColor());
	ratio->GetXaxis()->SetTitle("Detector-level Leading Jet p_{T}");
	ratio->GetYaxis()->SetTitle("Backfold/Measurement");
	ratio->Draw();

	TLine *l = new TLine(0,1,100,1);
	l->Draw();

		
}
